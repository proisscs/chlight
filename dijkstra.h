#ifndef FILE_DIJK_SEEN
#define FILE_DIJK_SEEN

#include <iostream>
#include <vector>
#include <random>
#include "c.h"
#include "Timer.h"
#include "CHParser.h"
#include <fstream>

using namespace std;

typedef pair<Rank, Distance> Edge;

struct Coordinates
{
    double lat;
    double lon;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar &lat;
        ar &lon;
    }
};

class Dijkstra
{
public:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar &_fwd_edges;
        ar &_bwd_edges;
        ar &_numEdges;
        ar &_osm_id_to_rank;
    }
    Dijkstra()
    {
    }

    void readFromFMIFile(string fname)
    {
        ifstream inputFile(fname.c_str());
        char junkC[1024];
        int junkI;
        string line;
        while (true)
        {
            getline(inputFile, line);
            if (line.find('#') == std::string::npos)
            {
                break;
            }
        }
        inputFile >> line;
        _numNodes = stoi(line);
        inputFile >> _numEdges;
        vector<pair<long, Rank>> osmID_rank_pairs(_numNodes);
        cout << "Num nodes: " << _numNodes << endl
             << "Num edges: " << _numEdges << endl;
        for (long i = 0; i < _numNodes; i++)
        {
            long long_tmp;
            double double_tmp;
            inputFile >> long_tmp >> long_tmp;
            osmID_rank_pairs[i] = make_pair(long_tmp, i);
            inputFile >> double_tmp >> double_tmp >> double_tmp;
        }
        _osm_id_to_rank = getOSMIDToRankMapping(osmID_rank_pairs);
        _fwd_edges.resize(_numNodes);
        _bwd_edges.resize(_numNodes);
        vector<EdgeType> edge_list(_numEdges);
        for (long j = 0; j < _numEdges; j++)
        {
            EdgeType new_edge;
            inputFile >> new_edge.source;
            inputFile >> new_edge.target;
            inputFile >> new_edge.weight;
            long tmp;
            inputFile >> tmp >> tmp;
            edge_list[j] = new_edge;
        }
        sort(edge_list.begin(), edge_list.end(), compareEdgeTypes);
        {
            if (_numEdges > 0)
            {
                const EdgeType &edge = edge_list[0];
                _fwd_edges[edge.source].push_back(make_pair(edge.target, edge.weight));
                _fwd_edges[edge.target].push_back(make_pair(edge.source, edge.weight));
            }
            for (long j = 1; j < edge_list.size(); j++)
            {
                const EdgeType &edge = edge_list[j];
                const EdgeType &prev_edge = edge_list[j - 1];
                if (edge.source != prev_edge.source || edge.target != prev_edge.target)
                {
                    _fwd_edges[edge.source].push_back(make_pair(edge.target, edge.weight));
                    _bwd_edges[edge.target].push_back(make_pair(edge.source, edge.weight));
                }
                else
                {
                    _numEdges--;
                }
            }
        }
        _fwdDistances.assign(_numNodes, c::NO_ENTRY);
        _bwdDistances.assign(_numNodes, c::NO_ENTRY);
        _added_by_bwd.resize(_numNodes);
        _added_by_fwd.resize(_numNodes);
        _is_initialized = true;
    }

    void constructDijkstraFromCH(string file)
    {
        vector<EdgeType> edge_list;
        vector<NodeID> rankToNodeID;
        Timer timer;
        {
            cout << "Loading CH..." << endl;
            CHParser ch_fmi;
            ch_fmi.readFromFMIFile(file);
            cout << "Finished." << endl;
            timer.start();
            _numNodes = ch_fmi.nofNodes();
            vector<Rank> nodeIDToRank;
            nodeIDToRank.resize(_numNodes);
            rankToNodeID = ch_fmi.rankToNodeID;
            vector<pair<long, Rank>> osm_rank_pairs(_numNodes);
            for (Rank i = 0; i < _numNodes; i++)
            {
                NodeID nodeID = rankToNodeID.at(i);
                osm_rank_pairs[i] = make_pair(ch_fmi.nodeList[nodeID].osmID, i);
            }
            _osm_id_to_rank = getOSMIDToRankMapping(osm_rank_pairs);
            for (Rank i = 0; i < _numNodes; i++)
            {
                Rank rank = _osm_id_to_rank[i];
                NodeID nodeID = rankToNodeID.at(rank);
                //DEBUG
                if(i==_numNodes/2)
                {
                    cout << "node ID for " << i << ": " << nodeID << endl;
                }
                //DEBUG END
                nodeIDToRank.at(nodeID) = i;
                _osm_id_to_rank[i] = i;
            }
            edge_list = ch_fmi.getEdges();
            for (EdgeType &edge : edge_list)
            {
                edge.source = nodeIDToRank[edge.source];
                edge.target = nodeIDToRank[edge.target];
            }
        }

        { // Save edges
            _fwd_edges.clear();
            _bwd_edges.clear();
            _fwd_edges.resize(_numNodes);
            _bwd_edges.resize(_numNodes);
            _numEdges = 0;

            sort(edge_list.begin(), edge_list.end(), compareEdgeTypes);
            cout << "Start initializing edges..." << endl;
            long shortcut_counter = 0;
            long current_rank = 0;
            for (const EdgeType &edge : edge_list)
            {
                if (edge.source == edge.target || edge.child_1 != c::NO_ENTRY)
                {
                    if (edge.child_1 != c::NO_ENTRY)
                    {
                        shortcut_counter++;
                    }
                }
                else
                {
                    _fwd_edges[edge.source].push_back(make_pair(edge.target, edge.weight));
                    _numEdges++;
                }
            }
            sort(edge_list.begin(), edge_list.end(), compareEdgeTypesBackwards);
            current_rank = 0;
            for (const EdgeType &edge : edge_list)
            {
                if (edge.target != edge.source && edge.child_1 == c::NO_ENTRY)
                {
                    _bwd_edges[edge.target].push_back(make_pair(edge.source, edge.weight));
                }
            }

            cout << "Construction finished.\nNumber of edges: " << _numEdges << "\nEdges CH: " << edge_list.size() << endl;
        }
        _fwdDistances.assign(_numNodes, c::NO_ENTRY);
        _bwdDistances.assign(_numNodes, c::NO_ENTRY);
        _added_by_bwd.resize(_numNodes);
        _added_by_fwd.resize(_numNodes);
        _is_initialized = true;
        timer.stop();
        cout << "Construction took " << timer.secs() << " seconds" << endl;
        //DEBUG
        for(int i=0; i<_fwd_edges[0].size(); i++)
        {
            cout << _fwd_edges[0][i].first << " " << _fwd_edges[0][i].second << endl;
        }
        //DEBUG END
    }

    vector<Edge> getShortestPathDijkstra(Rank source, Rank target)
    {
        typedef pair<Distance, Rank> PQElement;
        vector<Edge> shortest_path;
        if (source == target)
        {
            return shortest_path;
        }
        fill(_fwdDistances.begin(), _fwdDistances.end(), c::NO_ENTRY);
        fill(_bwdDistances.begin(), _bwdDistances.end(), c::NO_ENTRY);
        Distance best_dist = c::NO_ENTRY;
        Rank mid_point = c::NO_ENTRY;
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> fwd_pQ;
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> bwd_pQ;
        fwd_pQ.push(make_pair(0, source));
        bwd_pQ.push(make_pair(0, target));
        _fwdDistances[source] = 0;
        _bwdDistances[target] = 0;
        bool forward = true;
        bool should_extend = true;
        while (!fwd_pQ.empty() && !bwd_pQ.empty())
        {
            _counter++;
            if (forward)
            {
                PQElement next_element = fwd_pQ.top();
                Rank rank = next_element.second;
                Distance distance = next_element.first;
                fwd_pQ.pop();
                if (_fwdDistances[rank] == distance)
                {
                    if (_bwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = _bwdDistances[rank] + distance;
                        if (best_dist == c::NO_ENTRY || new_dist < best_dist)
                        {
                            best_dist = new_dist;
                            mid_point = rank;
                        }
                        if (_bwdDistances[rank] <= bwd_pQ.top().first)
                        {
                            should_extend = false;
                        }
                    }
                    if (should_extend)
                    {
                        for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                        {
                            const Edge &edge = _fwd_edges[rank][edgeid];
                            Distance new_distance = edge.second + distance;
                            if (_fwdDistances[edge.first] == c::NO_ENTRY || new_distance < _fwdDistances[edge.first])
                            {
                                _fwdDistances[edge.first] = new_distance;
                                _added_by_fwd[edge.first] = rank;
                                fwd_pQ.push(make_pair(new_distance, edge.first));
                            }
                        }
                    }
                }
            }
            else
            {
                PQElement next_element = bwd_pQ.top();
                Rank rank = next_element.second;
                Distance distance = next_element.first;
                bwd_pQ.pop();
                if (_bwdDistances[rank] == distance)
                {
                    if (_fwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = _fwdDistances[rank] + distance;
                        if (best_dist == c::NO_ENTRY || new_dist < best_dist)
                        {
                            best_dist = new_dist;
                            mid_point = rank;
                        }
                        if (_fwdDistances[rank] <= fwd_pQ.top().first)
                        {
                            should_extend = false;
                        }
                    }
                    if (should_extend)
                    {
                        for (int edgeid = 0; edgeid < _bwd_edges[rank].size(); edgeid++)
                        {
                            const Edge &edge = _bwd_edges[rank][edgeid];
                            Distance new_distance = edge.second + distance;
                            if (_bwdDistances[edge.first] == c::NO_ENTRY || new_distance < _bwdDistances[edge.first])
                            {
                                _bwdDistances[edge.first] = new_distance;
                                _added_by_bwd[edge.first] = rank;
                                bwd_pQ.push(make_pair(new_distance, edge.first));
                            }
                        }
                    }
                }
            }
            forward = !forward;
        }
        { // Construct shortest path
            if (mid_point != c::NO_ENTRY)
            {
                vector<Edge> fwd_shortest_path;
                vector<Edge> bwd_shortest_path;
                fwd_shortest_path.reserve(_numNodes / 1000);
                bwd_shortest_path.reserve(_numNodes / 1000);
                Rank cur_rank = mid_point;
                while (cur_rank != source)
                {
                    Rank prev_rank = _added_by_fwd[cur_rank];
                    for (const Edge &edge : _fwd_edges[prev_rank])
                    {
                        if (edge.first == cur_rank)
                        {
                            fwd_shortest_path.push_back(edge);
                            break;
                        }
                    }
                    cur_rank = prev_rank;
                }
                cur_rank = mid_point;
                while (cur_rank != target)
                {
                    Rank prev_rank = _added_by_bwd[cur_rank];
                    for (const Edge &edge : _bwd_edges[prev_rank])
                    {
                        if (edge.first == cur_rank)
                        {
                            bwd_shortest_path.push_back(make_pair(prev_rank, edge.second));
                            break;
                        }
                    }
                    cur_rank = prev_rank;
                }
                shortest_path.reserve(fwd_shortest_path.size() + bwd_shortest_path.size());
                for (int edgeid = fwd_shortest_path.size() - 1; edgeid >= 0; edgeid--)
                {
                    shortest_path.push_back(fwd_shortest_path[edgeid]);
                }
                shortest_path.insert(shortest_path.end(), bwd_shortest_path.begin(), bwd_shortest_path.end());
            }
        }
        return shortest_path;
    }

    Distance oneToAllDijkstra(Rank source)
    {
        // DEBUG
        Timer timer;
        timer.start();
        // DEBUG END
        typedef pair<Distance, Rank> PQElement;
        vector<Edge> shortest_path;
        fill(_fwdDistances.begin(), _fwdDistances.end(), c::NO_ENTRY);
        vector<PQElement> container;
        container.reserve(_numNodes);
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> pQ(greater<PQElement>(), move(container));
        pQ.push(make_pair(0, source));
        _fwdDistances[source] = 0;
        while (!pQ.empty())
        {
            _counter++;
            PQElement next_element = pQ.top();
            Rank rank = next_element.second;
            Distance distance = next_element.first;
            pQ.pop();
            if (_fwdDistances[rank] == distance)
            {
                for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                {
                    const Edge &edge = _fwd_edges[rank][edgeid];
                    Distance new_distance = edge.second + distance;
                    if (_fwdDistances[edge.first] == c::NO_ENTRY || new_distance < _fwdDistances[edge.first])
                    {
                        _fwdDistances[edge.first] = new_distance;
                        _added_by_fwd[edge.first] = rank;
                        pQ.push(make_pair(new_distance, edge.first));
                    }
                }
            }
        }
        // DEBUG
        timer.stop();
        cout << "Elapsed time: " << timer.secs() << endl;
        // DEBUG END
        return _fwdDistances[0];
    }

    long getNumberOfVisitedNodes(Rank source, Rank target)
    {
        long visited_nodes = 0;
        typedef pair<Distance, Rank> PQElement;
        vector<Edge> shortest_path;
        if (source == target)
        {
            return 1;
        }
        fill(_fwdDistances.begin(), _fwdDistances.end(), c::NO_ENTRY);
        fill(_bwdDistances.begin(), _bwdDistances.end(), c::NO_ENTRY);
        Distance best_dist = c::NO_ENTRY;
        Rank mid_point = c::NO_ENTRY;
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> fwd_pQ;
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> bwd_pQ;
        fwd_pQ.push(make_pair(0, source));
        bwd_pQ.push(make_pair(0, target));
        _fwdDistances[source] = 0;
        _bwdDistances[target] = 0;
        visited_nodes = 2;
        bool forward = true;
        while (!fwd_pQ.empty() && !bwd_pQ.empty())
        {
            _counter++;
            if (forward)
            {
                PQElement next_element = fwd_pQ.top();
                Rank rank = next_element.second;
                Distance distance = next_element.first;
                fwd_pQ.pop();
                if (_fwdDistances[rank] == distance)
                {
                    if (_bwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = _bwdDistances[rank] + distance;
                        if (best_dist == c::NO_ENTRY || new_dist < best_dist)
                        {
                            best_dist = new_dist;
                            mid_point = rank;
                        }
                    }
                    for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                    {
                        const Edge &edge = _fwd_edges[rank][edgeid];
                        Distance new_distance = edge.second + distance;
                        if (_fwdDistances[edge.first] == c::NO_ENTRY || new_distance < _fwdDistances[edge.first])
                        {
                            if (_fwdDistances[edge.first] == c::NO_ENTRY)
                            {
                                visited_nodes++;
                            }
                            _fwdDistances[edge.first] = new_distance;
                            if (best_dist == c::NO_ENTRY)
                            {
                                fwd_pQ.push(make_pair(new_distance, edge.first));
                            }
                        }
                    }
                }
            }
            else
            {
                PQElement next_element = bwd_pQ.top();
                Rank rank = next_element.second;
                Distance distance = next_element.first;
                bwd_pQ.pop();
                if (_bwdDistances[rank] == distance)
                {
                    if (_fwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = _fwdDistances[rank] + distance;
                        if (best_dist == c::NO_ENTRY || new_dist < best_dist)
                        {
                            best_dist = new_dist;
                            mid_point = rank;
                        }
                    }
                    for (int edgeid = 0; edgeid < _bwd_edges[rank].size(); edgeid++)
                    {
                        const Edge &edge = _bwd_edges[rank][edgeid];
                        Distance new_distance = edge.second + distance;
                        if (_bwdDistances[edge.first] == c::NO_ENTRY || new_distance < _bwdDistances[edge.first])
                        {
                            if (_bwdDistances[edge.first] == c::NO_ENTRY)
                            {
                                visited_nodes++;
                            }
                            _bwdDistances[edge.first] = new_distance;
                            if (best_dist == c::NO_ENTRY)
                            {
                                bwd_pQ.push(make_pair(new_distance, edge.first));
                            }
                        }
                    }
                }
            }
            forward = !forward;
        }
        return visited_nodes;
    }

    vector<Edge> getShortestPathUnidirectionalDijkstra(Rank source, Rank target)
    {
        typedef pair<Distance, Rank> PQElement;
        vector<Edge> shortest_path;
        if (source == target)
        {
            return shortest_path;
        }
        fill(_fwdDistances.begin(), _fwdDistances.end(), c::NO_ENTRY);
        // vector<bool> is_settled(_numNodes, false);
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> pQ;
        pQ.push(make_pair(0, source));
        _fwdDistances[source] = 0;
        while (!pQ.empty() && pQ.top().second != target)
        {
            _counter++;
            PQElement next_element = pQ.top();
            Rank rank = next_element.second;
            Distance distance = next_element.first;
            pQ.pop();
            if (_fwdDistances[rank] == distance)
            {
                for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                {
                    const Edge &edge = _fwd_edges[rank][edgeid];
                    Distance new_distance = edge.second + distance;
                    if (_fwdDistances[edge.first] == c::NO_ENTRY || new_distance < _fwdDistances[edge.first])
                    {
                        _fwdDistances[edge.first] = new_distance;
                        _added_by_fwd[edge.first] = rank;
                        pQ.push(make_pair(new_distance, edge.first));
                    }
                }
            }
        }
        { // Construct shortest path
            Distance path_dist = 0;
            if (!pQ.empty())
            {
                vector<Edge> fwd_shortest_path;
                fwd_shortest_path.reserve(_numNodes / 100);
                Rank cur_rank = target;
                while (cur_rank != source)
                {
                    Rank prev_rank = _added_by_fwd[cur_rank];
                    for (const Edge &edge : _fwd_edges[prev_rank])
                    {
                        if (edge.first == cur_rank)
                        {
                            fwd_shortest_path.push_back(edge);
                            path_dist += edge.second;
                            break;
                        }
                    }
                    cur_rank = prev_rank;
                }
                shortest_path.reserve(fwd_shortest_path.size());
                for (int edgeid = fwd_shortest_path.size() - 1; edgeid >= 0; edgeid--)
                {
                    shortest_path.push_back(fwd_shortest_path[edgeid]);
                }
            }
        }

        return shortest_path;
    }

    Distance getDistanceUnidirectionalDijkstra(Rank source, Rank target)
    {
        typedef pair<Distance, Rank> PQElement;
        if (source == target)
        {
            return 0;
        }
        fill(_fwdDistances.begin(), _fwdDistances.end(), c::NO_ENTRY);
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> pQ;
        pQ.push(make_pair(0, source));
        _fwdDistances[source] = 0;
        while (!pQ.empty() && pQ.top().second != target)
        {
            PQElement next_element = pQ.top();
            Rank rank = next_element.second;
            Distance distance = next_element.first;
            pQ.pop();
            if (_fwdDistances[rank] == distance)
            {
                for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                {
                    const Edge &edge = _fwd_edges[rank][edgeid];
                    Distance new_distance = edge.second + distance;
                    if (_fwdDistances[edge.first] == c::NO_ENTRY || new_distance < _fwdDistances[edge.first])
                    {
                        _fwdDistances[edge.first] = new_distance;
                        pQ.push(make_pair(new_distance, edge.first));
                    }
                }
            }
        }
        return _fwdDistances[target];
    }

    vector<Rank> getNodesOfDijkstraRanks(Rank source, const vector<long> &dijksta_ranks)
    {
        typedef pair<Distance, Rank> PQElement;
        fill(_fwdDistances.begin(), _fwdDistances.end(), c::NO_ENTRY);
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> pQ;
        pQ.push(make_pair(0, source));
        _fwdDistances[source] = 0;
        int settled_nodes_counter = 0;
        int saved_nodes_counter = 0;
        vector<Rank> settled_nodes(dijksta_ranks.size(), c::NO_ENTRY);
        if (saved_nodes_counter < dijksta_ranks.size() && dijksta_ranks[saved_nodes_counter] == settled_nodes_counter)
        {
            settled_nodes[saved_nodes_counter++] = source;
        }
        while (!pQ.empty() && saved_nodes_counter < dijksta_ranks.size())
        {
            PQElement next_element = pQ.top();
            Rank rank = next_element.second;
            Distance distance = next_element.first;
            pQ.pop();
            if (_fwdDistances[rank] == distance)
            {
                settled_nodes_counter++;
                if (saved_nodes_counter < dijksta_ranks.size() && dijksta_ranks[saved_nodes_counter] == settled_nodes_counter)
                {
                    settled_nodes[saved_nodes_counter++] = rank;
                }
                for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                {
                    const Edge &edge = _fwd_edges[rank][edgeid];
                    Distance new_distance = edge.second + distance;
                    if (_fwdDistances[edge.first] == c::NO_ENTRY || new_distance < _fwdDistances[edge.first])
                    {
                        _fwdDistances[edge.first] = new_distance;
                        pQ.push(make_pair(new_distance, edge.first));
                    }
                }
            }
        }
        return settled_nodes;
    }

    long nofNodes()
    {
        return _numNodes;
    }

    long nofEdges()
    {
        return _numEdges;
    }

    long getSizeOfDijkstra()
    {
        return 2 * _numEdges * sizeof(Edge) + 2 * sizeof(Rank) * _numNodes;
    }

    long getSizeOfUniDijkstra()
    {
        return _numEdges * sizeof(Edge) + sizeof(Rank) * _numNodes;
    }

    long getSizeOfAStar()
    {
        return _numEdges * sizeof(Edge) + sizeof(Rank) * _numNodes + sizeof(Coordinates) * _numNodes;
    }

    void resetCounter()
    {
        _counter = 0;
    }

    long getCounter()
    {
        return _counter;
    }

    bool isInitialized()
    {
        return _is_initialized;
    }

    Rank osmIDToRank(long osmID)
    {
        return _osm_id_to_rank[osmID];
    }

private:
    long _numNodes;
    long _numEdges;
    long _counter = 0;
    bool _is_initialized = false;
    vector<vector<Edge>> _fwd_edges;
    vector<vector<Edge>> _bwd_edges;

    vector<Distance> _fwdDistances;
    vector<Distance> _bwdDistances;
    vector<Rank> _added_by_fwd;
    vector<Rank> _added_by_bwd;
    vector<Rank> _osm_id_to_rank;
};

#endif