#ifndef FILE_CH_SEEN
#define FILE_CH_SEEN

#include <iostream>
#include <vector>
#include <random>
#include "c.h"
#include "Timer.h"
#include "CHParser.h"
#include <fstream>

using namespace std;

class CH
{
public:
    CH()
    {
    }
    void constructCH(string file)
    {
        vector<EdgeType> edge_list;
        vector<NodeID> rankToNodeID;
        vector<Rank> nodeIDToRank;
        {
            cout << "Loading CH..." << endl;
            CHParser ch_fmi;
            ch_fmi.readFromFMIFile(file);
            cout << "Finished." << endl;
            _numNodes = ch_fmi.nofNodes();
            nodeIDToRank.resize(_numNodes);
            rankToNodeID = ch_fmi.rankToNodeID;
            vector<pair<long,Rank>> osm_rank_pairs(_numNodes);
            for (Rank i = 0; i < _numNodes; i++)
            {
                NodeID nodeID = rankToNodeID.at(i);
                nodeIDToRank.at(nodeID) = i;
                osm_rank_pairs[i] = make_pair(ch_fmi.nodeList[nodeID].osmID, i);
            }
            _osm_id_to_rank = getOSMIDToRankMapping(osm_rank_pairs);
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
            _fwd_edges.reserve(edge_list.size());
            _bwd_edges.reserve(edge_list.size());
            _numEdges = 0;
            _numShortcuts = 0;
            int count_unnecessary_edges = 0;
            cout << "Start initializing edges..." << endl;
            sort(edge_list.begin(), edge_list.end(), compareEdgeTypes);
            vector<EdgeID> index_to_edgeid(edge_list.size());
            {
                long counter = 0;
                for (const EdgeType &edge : edge_list)
                {
                    index_to_edgeid[edge.index] = counter;
                    counter++;
                }
            }
            _fwd_shortcuts.resize(_numNodes);
            _fwd_edges.resize(_numNodes);
            for (EdgeID edgeid = 0; edgeid < edge_list.size(); edgeid++)
            {
                const EdgeType &edge = edge_list[edgeid];
                if (edge.source < edge.target)
                {
                    if (edgeid == 0 || edge_list[edgeid - 1].source != edge.source || edge_list[edgeid - 1].target != edge.target)
                    {
                        if (edge.child_1 == c::NO_ENTRY)
                        {
                            Edge new_edge = make_pair(edge.target, edge.weight);
                            _fwd_edges.at(edge.source).push_back(new_edge);
                            _numEdges++;
                        }
                        else
                        {
                            Shortcut new_shortcut = {
                                edge.target,
                                edge.weight,
                                edge_list.at(index_to_edgeid[edge.child_1]).target,
                            };
                            _fwd_shortcuts.at(edge.source).push_back(new_shortcut);
                            _numShortcuts++;
                        }
                    }
                    else
                    {
                        count_unnecessary_edges++;
                    }
                }
            }
            sort(edge_list.begin(), edge_list.end(), compareEdgeTypesBackwards);
            {
                long counter = 0;
                for (const EdgeType &edge : edge_list)
                {
                    index_to_edgeid[edge.index] = counter;
                    counter++;
                }
            }
            _bwd_shortcuts.resize(_numNodes);
            _bwd_edges.resize(_numNodes);
            for (EdgeID edgeid = 0; edgeid < edge_list.size(); edgeid++)
            {
                const EdgeType &edge = edge_list[edgeid];
                if (edge.source > edge.target)
                {
                    if (edgeid == 0 || edge_list[edgeid - 1].source != edge.source || edge_list[edgeid - 1].target != edge.target)
                    {
                        if (edge.child_1 == c::NO_ENTRY)
                        {
                            Edge new_edge = make_pair(edge.source, edge.weight);
                            _bwd_edges.at(edge.target).push_back(new_edge);
                            _numEdges++;
                        }
                        else
                        {
                            Shortcut new_shortcut = {
                                edge.source,
                                edge.weight,
                                edge_list.at(index_to_edgeid[edge.child_1]).target,
                            };
                            _bwd_shortcuts.at(edge.target).push_back(new_shortcut);
                            _numShortcuts++;
                        }
                    }
                    else
                    {
                        count_unnecessary_edges++;
                    }
                }
            }
            cout << "Construction finished.\nNumber of edges: " << _numEdges << "\nEdges original CH: " << edge_list.size() << "\nUnnecessary edges: " << count_unnecessary_edges << endl;
        }
        _fwdDistances.assign(_numNodes, c::NO_ENTRY);
        _bwdDistances.assign(_numNodes, c::NO_ENTRY);
        _added_by_bwd.resize(_numNodes);
        _added_by_fwd.resize(_numNodes);
        _is_initialized = true;
    }

    bool isInitialized()
    {
        return _is_initialized;
    }

    Distance getDistance(Rank source, Rank target)
    {
        typedef pair<Distance, Rank> PQElement;
        if (source == target)
        {
            return 0;
        }
        vector<Rank> fwd_visited_nodes;
        fwd_visited_nodes.reserve(10000);
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> pQ;
        pQ.push(make_pair(0, source));
        fwd_visited_nodes.push_back(source);
        _fwdDistances.at(source) = 0;
        while (!pQ.empty())
        {
            PQElement next_element = pQ.top();
            Rank rank = next_element.second;
            Distance distance = next_element.first;
            pQ.pop();
            if (_fwdDistances[rank] == distance)
            {
                bool should_be_stalled = false;
                for (const Shortcut &edge : _bwd_shortcuts[rank])
                {
                    if (_fwdDistances[edge.node] != c::NO_ENTRY && _fwdDistances[edge.node] + edge.distance < distance)
                    {
                        should_be_stalled = true;
                        break;
                    }
                }
                if (!should_be_stalled)
                {
                    for (const Edge &edge : _bwd_edges[rank])
                    {
                        if (_fwdDistances[edge.first] != c::NO_ENTRY && _fwdDistances[edge.first] + edge.second < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                }
                if (!should_be_stalled)
                {
                    for (const Shortcut &edge : _fwd_shortcuts[rank])
                    {
                        Distance new_dist = distance + edge.distance;
                        if (_fwdDistances[edge.node] == c::NO_ENTRY || _fwdDistances[edge.node] > new_dist)
                        {
                            if (_fwdDistances[edge.node] == c::NO_ENTRY)
                            {
                                fwd_visited_nodes.push_back(edge.node);
                            }
                            _fwdDistances[edge.node] = new_dist;
                            pQ.push(make_pair(new_dist, edge.node));
                        }
                    }
                    for (const Edge &edge : _fwd_edges[rank])
                    {
                        Distance new_dist = distance + edge.second;
                        if (_fwdDistances[edge.first] == c::NO_ENTRY || _fwdDistances[edge.first] > new_dist)
                        {
                            if (_fwdDistances[edge.first] == c::NO_ENTRY)
                            {
                                fwd_visited_nodes.push_back(edge.first);
                            }
                            _fwdDistances[edge.first] = new_dist;
                            pQ.push(make_pair(new_dist, edge.first));
                        }
                    }
                }
            }
        }
        pQ.push(make_pair(0, target));
        _bwdDistances[target] = 0;
        vector<Rank> bwd_visited_nodes;
        bwd_visited_nodes.reserve(10000);
        bwd_visited_nodes.push_back(target);
        while (!pQ.empty())
        {
            PQElement next_element = pQ.top();
            Rank rank = next_element.second;
            Distance distance = next_element.first;
            pQ.pop();
            if (_bwdDistances[rank] == distance)
            {
                bool should_be_stalled = false;
                for (const Shortcut &edge : _fwd_shortcuts[rank])
                {
                    if (_bwdDistances[edge.node] != c::NO_ENTRY && _bwdDistances[edge.node] + edge.distance < distance)
                    {
                        should_be_stalled = true;
                        break;
                    }
                }
                if (!should_be_stalled)
                {
                    for (const Edge &edge : _fwd_edges[rank])
                    {
                        if (_bwdDistances[edge.first] != c::NO_ENTRY && _bwdDistances[edge.first] + edge.second < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                }
                if (!should_be_stalled)
                {
                    for (const Shortcut &edge : _bwd_shortcuts[rank])
                    {
                        Distance new_dist = distance + edge.distance;
                        if (_bwdDistances[edge.node] == c::NO_ENTRY || _bwdDistances[edge.node] > new_dist)
                        {
                            if (_bwdDistances[edge.node] == c::NO_ENTRY)
                            {
                                bwd_visited_nodes.push_back(edge.node);
                            }
                            _bwdDistances[edge.node] = new_dist;
                            pQ.push(make_pair(new_dist, edge.node));
                        }
                    }
                    for (const Edge &edge : _bwd_edges[rank])
                    {
                        Distance new_dist = distance + edge.second;
                        if (_bwdDistances[edge.first] == c::NO_ENTRY || _bwdDistances[edge.first] > new_dist)
                        {
                            if (_bwdDistances[edge.first] == c::NO_ENTRY)
                            {
                                bwd_visited_nodes.push_back(edge.first);
                            }
                            _bwdDistances[edge.first] = new_dist;
                            pQ.push(make_pair(new_dist, edge.first));
                        }
                    }
                }
            }
        }
        Distance best_dist = c::NO_ENTRY;
        for (auto rank : fwd_visited_nodes)
        {
            if (_bwdDistances[rank] != c::NO_ENTRY)
            {
                Distance new_dist = _fwdDistances[rank] + _bwdDistances[rank];
                if (best_dist == c::NO_ENTRY || best_dist > new_dist)
                {
                    best_dist = new_dist;
                }
            }
            _fwdDistances[rank] = c::NO_ENTRY;
        }
        for (auto rank : bwd_visited_nodes)
        {
            _bwdDistances[rank] = c::NO_ENTRY;
        }
        return best_dist;
    }

    long getNumberVisitedNodes(Rank source, Rank target)
    {   
        long num_visited = 0;
        typedef pair<Distance, Rank> PQElement;
        if (source == target)
        {
            return 0;
        }
        vector<Rank> fwd_visited_nodes;
        fwd_visited_nodes.reserve(10000);
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> pQ;
        pQ.push(make_pair(0, source));
        fwd_visited_nodes.push_back(source);
        _fwdDistances.at(source) = 0;
        while (!pQ.empty())
        {
            PQElement next_element = pQ.top();
            Rank rank = next_element.second;
            Distance distance = next_element.first;
            pQ.pop();
            if (_fwdDistances[rank] == distance)
            {
                bool should_be_stalled = false;
                for (const Shortcut &edge : _bwd_shortcuts[rank])
                {
                    if (_fwdDistances[edge.node] != c::NO_ENTRY && _fwdDistances[edge.node] + edge.distance < distance)
                    {
                        should_be_stalled = true;
                        break;
                    }
                }
                if (!should_be_stalled)
                {
                    for (const Edge &edge : _bwd_edges[rank])
                    {
                        if (_fwdDistances[edge.first] != c::NO_ENTRY && _fwdDistances[edge.first] + edge.second < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                }
                if (!should_be_stalled)
                {
                    for (const Shortcut &edge : _fwd_shortcuts[rank])
                    {
                        Distance new_dist = distance + edge.distance;
                        if (_fwdDistances[edge.node] == c::NO_ENTRY || _fwdDistances[edge.node] > new_dist)
                        {
                            if (_fwdDistances[edge.node] == c::NO_ENTRY)
                            {
                                fwd_visited_nodes.push_back(edge.node);
                            }
                            _fwdDistances[edge.node] = new_dist;
                            pQ.push(make_pair(new_dist, edge.node));
                        }
                    }
                    for (const Edge &edge : _fwd_edges[rank])
                    {
                        Distance new_dist = distance + edge.second;
                        if (_fwdDistances[edge.first] == c::NO_ENTRY || _fwdDistances[edge.first] > new_dist)
                        {
                            if (_fwdDistances[edge.first] == c::NO_ENTRY)
                            {
                                fwd_visited_nodes.push_back(edge.first);
                            }
                            _fwdDistances[edge.first] = new_dist;
                            pQ.push(make_pair(new_dist, edge.first));
                        }
                    }
                }
            }
        }
        pQ.push(make_pair(0, target));
        _bwdDistances[target] = 0;
        vector<Rank> bwd_visited_nodes;
        bwd_visited_nodes.reserve(10000);
        bwd_visited_nodes.push_back(target);
        while (!pQ.empty())
        {
            PQElement next_element = pQ.top();
            Rank rank = next_element.second;
            Distance distance = next_element.first;
            pQ.pop();
            if (_bwdDistances[rank] == distance)
            {
                bool should_be_stalled = false;
                for (const Shortcut &edge : _fwd_shortcuts[rank])
                {
                    if (_bwdDistances[edge.node] != c::NO_ENTRY && _bwdDistances[edge.node] + edge.distance < distance)
                    {
                        should_be_stalled = true;
                        break;
                    }
                }
                if (!should_be_stalled)
                {
                    for (const Edge &edge : _fwd_edges[rank])
                    {
                        if (_bwdDistances[edge.first] != c::NO_ENTRY && _bwdDistances[edge.first] + edge.second < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                }
                if (!should_be_stalled)
                {
                    for (const Shortcut &edge : _bwd_shortcuts[rank])
                    {
                        Distance new_dist = distance + edge.distance;
                        if (_bwdDistances[edge.node] == c::NO_ENTRY || _bwdDistances[edge.node] > new_dist)
                        {
                            if (_bwdDistances[edge.node] == c::NO_ENTRY)
                            {
                                bwd_visited_nodes.push_back(edge.node);
                            }
                            _bwdDistances[edge.node] = new_dist;
                            pQ.push(make_pair(new_dist, edge.node));
                        }
                    }
                    for (const Edge &edge : _bwd_edges[rank])
                    {
                        Distance new_dist = distance + edge.second;
                        if (_bwdDistances[edge.first] == c::NO_ENTRY || _bwdDistances[edge.first] > new_dist)
                        {
                            if (_bwdDistances[edge.first] == c::NO_ENTRY)
                            {
                                bwd_visited_nodes.push_back(edge.first);
                            }
                            _bwdDistances[edge.first] = new_dist;
                            pQ.push(make_pair(new_dist, edge.first));
                        }
                    }
                }
            }
        }
        Distance best_dist = c::NO_ENTRY;
        for (auto rank : fwd_visited_nodes)
        {
            if (_bwdDistances[rank] != c::NO_ENTRY)
            {
                Distance new_dist = _fwdDistances[rank] + _bwdDistances[rank];
                if (best_dist == c::NO_ENTRY || best_dist > new_dist)
                {
                    best_dist = new_dist;
                }
            }
            _fwdDistances[rank] = c::NO_ENTRY;
        }
        for (auto rank : bwd_visited_nodes)
        {
            _bwdDistances[rank] = c::NO_ENTRY;
        }
        return fwd_visited_nodes.size() + bwd_visited_nodes.size();
    }

    pair<long, long> getNumNodesInUnpackedSearchTreeAndExpandedNodes(Rank source, Rank target)
    {
        typedef pair<Distance, Rank> PQElement;
        if (source == target)
        {
            return make_pair(0, 0);
        }
        vector<Rank> fwd_visited_nodes;
        fwd_visited_nodes.reserve(10000);
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> fwd_pQ;
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> bwd_pQ;
        fwd_pQ.push(make_pair(0, source));
        bwd_pQ.push(make_pair(0, target));
        fwd_visited_nodes.push_back(source);
        _fwdDistances[source] = 0;
        _bwdDistances[target] = 0;
        vector<Rank> bwd_visited_nodes;
        bwd_visited_nodes.reserve(10000);
        bwd_visited_nodes.push_back(target);
        Distance best_dist = c::NO_ENTRY;
        Rank midpoint = c::NO_ENTRY;
        struct TmpEdge {
            Rank source;
            Rank target;
            bool forward;
        };
        vector<TmpEdge> expanded_shortcuts;
        vector<TmpEdge> search_tree_shortcuts;
        bool forward = true;
        while (!fwd_pQ.empty() || !bwd_pQ.empty())
        {
            _counter++;
            if (forward)
            {
                PQElement next_element = fwd_pQ.top();
                Rank rank = next_element.second;
                Distance distance = next_element.first;
                if (best_dist != c::NO_ENTRY && distance > best_dist)
                {
                    break;
                }
                fwd_pQ.pop();
                if (_fwdDistances[rank] == distance)
                {
                    if (_bwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = distance + _bwdDistances[rank];
                        if (best_dist == c::NO_ENTRY || best_dist > new_dist)
                        {
                            best_dist = new_dist;
                            midpoint = rank;
                        }
                    }
                    bool should_be_stalled = false;
                    for (const Shortcut &edge : _bwd_shortcuts[rank])
                    {
                        if (_fwdDistances[edge.node] != c::NO_ENTRY && _fwdDistances[edge.node] + edge.distance < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                    if (!should_be_stalled)
                    {
                        for (const Edge &edge : _bwd_edges[rank])
                        {
                            if (_fwdDistances[edge.first] != c::NO_ENTRY && _fwdDistances[edge.first] + edge.second < distance)
                            {
                                should_be_stalled = true;
                                break;
                            }
                        }
                    }
                    if (!should_be_stalled)
                    {
                        for (int edgeid = 0; edgeid < _fwd_shortcuts[rank].size(); edgeid++)
                        {
                            const Shortcut &edge = _fwd_shortcuts[rank][edgeid];
                            Distance new_dist = distance + edge.distance;
                            if (_fwdDistances[edge.node] == c::NO_ENTRY || _fwdDistances[edge.node] > new_dist)
                            {
                                if (_fwdDistances[edge.node] == c::NO_ENTRY)
                                {
                                    fwd_visited_nodes.push_back(edge.node);
                                }
                                _fwdDistances[edge.node] = new_dist;
                                fwd_pQ.push(make_pair(new_dist, edge.node));
                                _added_by_fwd[edge.node] = rank;
                                expanded_shortcuts.push_back(TmpEdge {
                                    rank,
                                    edge.node,
                                    true
                                });
                            }
                        }
                        for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                        {
                            const Edge &edge = _fwd_edges[rank][edgeid];
                            Distance new_dist = distance + edge.second;
                            if (_fwdDistances[edge.first] == c::NO_ENTRY || _fwdDistances[edge.first] > new_dist)
                            {
                                if (_fwdDistances[edge.first] == c::NO_ENTRY)
                                {
                                    fwd_visited_nodes.push_back(edge.first);
                                }
                                _fwdDistances[edge.first] = new_dist;
                                fwd_pQ.push(make_pair(new_dist, edge.first));
                                _added_by_fwd[edge.first] = rank;
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
                if (best_dist != c::NO_ENTRY && distance > best_dist)
                {
                    break;
                }
                bwd_pQ.pop();
                if (_bwdDistances[rank] == distance)
                {
                    if (_fwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = distance + _fwdDistances[rank];
                        if (best_dist == c::NO_ENTRY || best_dist > new_dist)
                        {
                            best_dist = new_dist;
                            midpoint = rank;
                        }
                    }
                    bool should_be_stalled = false;
                    for (const Shortcut &edge : _fwd_shortcuts[rank])
                    {
                        if (_bwdDistances[edge.node] != c::NO_ENTRY && _bwdDistances[edge.node] + edge.distance < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                    if (!should_be_stalled)
                    {
                        for (const Edge &edge : _fwd_edges[rank])
                        {
                            if (_bwdDistances[edge.first] != c::NO_ENTRY && _bwdDistances[edge.first] + edge.second < distance)
                            {
                                should_be_stalled = true;
                                break;
                            }
                        }
                    }
                    if (!should_be_stalled)
                    {
                        for (int edgeid = 0; edgeid < _bwd_shortcuts[rank].size(); edgeid++)
                        {
                            const Shortcut &edge = _bwd_shortcuts[rank][edgeid];
                            Distance new_dist = distance + edge.distance;
                            if (_bwdDistances[edge.node] == c::NO_ENTRY || _bwdDistances[edge.node] > new_dist)
                            {
                                if (_bwdDistances[edge.node] == c::NO_ENTRY)
                                {
                                    bwd_visited_nodes.push_back(edge.node);
                                }
                                _bwdDistances[edge.node] = new_dist;
                                bwd_pQ.push(make_pair(new_dist, edge.node));
                                _added_by_bwd[edge.node] = rank;
                                expanded_shortcuts.push_back(TmpEdge {
                                    edge.node,
                                    rank,
                                    false
                                });
                            }
                        }
                        for (int edgeid = 0; edgeid < _bwd_edges[rank].size(); edgeid++)
                        {
                            const Edge &edge = _bwd_edges[rank][edgeid];
                            Distance new_dist = distance + edge.second;
                            if (_bwdDistances[edge.first] == c::NO_ENTRY || _bwdDistances[edge.first] > new_dist)
                            {
                                if (_bwdDistances[edge.first] == c::NO_ENTRY)
                                {
                                    bwd_visited_nodes.push_back(edge.first);
                                }
                                _bwdDistances[edge.first] = new_dist;
                                bwd_pQ.push(make_pair(new_dist, edge.first));
                                _added_by_bwd[edge.first] = rank;
                            }
                        }
                    }
                }
            }
            if (fwd_pQ.empty())
            {
                forward = false;
            }
            else if (bwd_pQ.empty())
            {
                forward = true;
            }
            else
            {
                forward = fwd_pQ.top().first < bwd_pQ.top().first;
            }
        }
        for (auto rank : fwd_visited_nodes)
        {
            _fwdDistances[rank] = c::NO_ENTRY;
            if (rank != source)
            {
                search_tree_shortcuts.push_back(TmpEdge {
                                    _added_by_fwd[rank],
                                    rank,
                                    true
                                });
            }
        }
        for (auto rank : bwd_visited_nodes)
        {
            _bwdDistances[rank] = c::NO_ENTRY;
            if (rank != target)
            {
                search_tree_shortcuts.push_back(TmpEdge {
                                    rank,
                                    _added_by_bwd[rank],
                                    false
                                });
            }
        }
        long count_visited_tree = 0;
        vector<bool> is_visited_fwd(_numNodes, false);
        vector<bool> is_visited_bwd(_numNodes, false);
        for (int i = 0; i < search_tree_shortcuts.size(); i++)
        {
            Rank edge_source = search_tree_shortcuts[i].source;
            Rank edge_target = search_tree_shortcuts[i].target;
            bool edge_forward = search_tree_shortcuts[i].forward;
            bool is_forward = edge_source < edge_target;
            if(edge_forward)
            {
                if(!is_visited_fwd[edge_source])
                {
                    is_visited_fwd[edge_source] = true;
                    count_visited_tree++;
                }
                if(!is_visited_fwd[edge_target])
                {
                    is_visited_fwd[edge_target] = true;
                    count_visited_tree++;
                }
            }
            else
            {
                if(!is_visited_bwd[edge_source])
                {
                    is_visited_bwd[edge_source] = true;
                    count_visited_tree++;
                }
                if(!is_visited_bwd[edge_target])
                {
                    is_visited_bwd[edge_target] = true;
                    count_visited_tree++;
                }
            }
            if (is_forward)
            {
                for (const Shortcut &edge : _fwd_shortcuts[edge_source])
                {
                    if (edge.node == edge_target)
                    {
                        Rank edge_midnode = edge.midnode;
                        search_tree_shortcuts.push_back(TmpEdge {
                            edge_source,
                            edge_midnode,
                            edge_forward
                        });
                        search_tree_shortcuts.push_back(TmpEdge {
                            edge_midnode,
                            edge_target,
                            edge_forward
                        });
                        break;
                    }
                }
            }
            else
            {
                for (const Shortcut &edge : _bwd_shortcuts[edge_target])
                {
                    if (edge.node == edge_source)
                    {
                        Rank edge_midnode = edge.midnode;
                        search_tree_shortcuts.push_back(TmpEdge {
                            edge_source,
                            edge_midnode,
                            edge_forward
                        });
                        search_tree_shortcuts.push_back(TmpEdge {
                            edge_midnode,
                            edge_target,
                            edge_forward
                        });
                        break;
                    }
                }
            }
        }
        long count_visited_expanded = 0;
        fill(is_visited_fwd.begin(), is_visited_fwd.end(), false);
        fill(is_visited_bwd.begin(), is_visited_bwd.end(), false);
        for (int i = 0; i < expanded_shortcuts.size(); i++)
        {
            Rank edge_source = expanded_shortcuts[i].source;
            Rank edge_target = expanded_shortcuts[i].target;
            bool edge_forward = expanded_shortcuts[i].forward;
            if(edge_forward)
            {
                if(!is_visited_fwd[edge_source])
                {
                    is_visited_fwd[edge_source] = true;
                    count_visited_expanded++;
                }
                if(!is_visited_fwd[edge_target])
                {
                    is_visited_fwd[edge_target] = true;
                    count_visited_expanded++;
                }
            }
            else
            {
                if(!is_visited_bwd[edge_source])
                {
                    is_visited_bwd[edge_source] = true;
                    count_visited_expanded++;
                }
                if(!is_visited_bwd[edge_target])
                {
                    is_visited_bwd[edge_target] = true;
                    count_visited_expanded++;
                }
            }
            bool is_forward = edge_source < edge_target;
            if (is_forward)
            {
                for (const Shortcut &edge : _fwd_shortcuts[edge_source])
                {
                    if (edge.node == edge_target)
                    {
                        Rank edge_midnode = edge.midnode;
                        expanded_shortcuts.push_back(TmpEdge {
                            edge_source,
                            edge_midnode,
                            edge_forward
                        });
                        expanded_shortcuts.push_back(TmpEdge {
                            edge_midnode,
                            edge_target,
                            edge_forward
                        });
                        break;
                    }
                }
            }
            else
            {
                for (const Shortcut &edge : _bwd_shortcuts[edge_target])
                {
                    if (edge.node == edge_source)
                    {
                        Rank edge_midnode = edge.midnode;
                        expanded_shortcuts.push_back(TmpEdge {
                            edge_source,
                            edge_midnode,
                            edge_forward
                        });
                        expanded_shortcuts.push_back(TmpEdge {
                            edge_midnode,
                            edge_target,
                            edge_forward
                        });
                        break;
                    }
                }
            }
        }
        return make_pair(count_visited_tree, count_visited_expanded);
    }

    vector<Edge> getShortestPath(Rank source, Rank target)
    {
        typedef pair<Distance, Rank> PQElement;
        vector<Edge> shortest_path;
        if (source == target)
        {
            return shortest_path;
        }
        vector<Rank> fwd_visited_nodes;
        fwd_visited_nodes.reserve(10000);
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> fwd_pQ;
        priority_queue<PQElement, vector<PQElement>, greater<PQElement>> bwd_pQ;
        fwd_pQ.push(make_pair(0, source));
        bwd_pQ.push(make_pair(0, target));
        fwd_visited_nodes.push_back(source);
        _fwdDistances[source] = 0;
        _bwdDistances[target] = 0;
        vector<Rank> bwd_visited_nodes;
        bwd_visited_nodes.reserve(10000);
        bwd_visited_nodes.push_back(target);
        Distance best_dist = c::NO_ENTRY;
        Rank midpoint = c::NO_ENTRY;
        bool forward = true;
        while (!fwd_pQ.empty() || !bwd_pQ.empty())
        {
            _counter++;
            if (forward)
            {
                PQElement next_element = fwd_pQ.top();
                Rank rank = next_element.second;
                Distance distance = next_element.first;
                if (best_dist != c::NO_ENTRY && distance > best_dist)
                {
                    break;
                }
                fwd_pQ.pop();
                if (_fwdDistances[rank] == distance)
                {
                    if (_bwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = distance + _bwdDistances[rank];
                        if (best_dist == c::NO_ENTRY || best_dist > new_dist)
                        {
                            best_dist = new_dist;
                            midpoint = rank;
                        }
                    }
                    bool should_be_stalled = false;
                    for (const Shortcut &edge : _bwd_shortcuts[rank])
                    {
                        if (_fwdDistances[edge.node] != c::NO_ENTRY && _fwdDistances[edge.node] + edge.distance < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                    if (!should_be_stalled)
                    {
                        for (const Edge &edge : _bwd_edges[rank])
                        {
                            if (_fwdDistances[edge.first] != c::NO_ENTRY && _fwdDistances[edge.first] + edge.second < distance)
                            {
                                should_be_stalled = true;
                                break;
                            }
                        }
                    }
                    if (!should_be_stalled)
                    {
                        for (int edgeid = 0; edgeid < _fwd_shortcuts[rank].size(); edgeid++)
                        {
                            const Shortcut &edge = _fwd_shortcuts[rank][edgeid];
                            Distance new_dist = distance + edge.distance;
                            if (_fwdDistances[edge.node] == c::NO_ENTRY || _fwdDistances[edge.node] > new_dist)
                            {
                                if (_fwdDistances[edge.node] == c::NO_ENTRY)
                                {
                                    fwd_visited_nodes.push_back(edge.node);
                                }
                                _fwdDistances[edge.node] = new_dist;
                                fwd_pQ.push(make_pair(new_dist, edge.node));
                                _added_by_fwd[edge.node] = rank;
                            }
                        }
                        for (int edgeid = 0; edgeid < _fwd_edges[rank].size(); edgeid++)
                        {
                            const Edge &edge = _fwd_edges[rank][edgeid];
                            Distance new_dist = distance + edge.second;
                            if (_fwdDistances[edge.first] == c::NO_ENTRY || _fwdDistances[edge.first] > new_dist)
                            {
                                if (_fwdDistances[edge.first] == c::NO_ENTRY)
                                {
                                    fwd_visited_nodes.push_back(edge.first);
                                }
                                _fwdDistances[edge.first] = new_dist;
                                fwd_pQ.push(make_pair(new_dist, edge.first));
                                _added_by_fwd[edge.first] = rank;
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
                if (best_dist != c::NO_ENTRY && distance > best_dist)
                {
                    break;
                }
                bwd_pQ.pop();
                if (_bwdDistances[rank] == distance)
                {
                    if (_fwdDistances[rank] != c::NO_ENTRY)
                    {
                        Distance new_dist = distance + _fwdDistances[rank];
                        if (best_dist == c::NO_ENTRY || best_dist > new_dist)
                        {
                            best_dist = new_dist;
                            midpoint = rank;
                        }
                    }
                    bool should_be_stalled = false;
                    for (const Shortcut &edge : _fwd_shortcuts[rank])
                    {
                        if (_bwdDistances[edge.node] != c::NO_ENTRY && _bwdDistances[edge.node] + edge.distance < distance)
                        {
                            should_be_stalled = true;
                            break;
                        }
                    }
                    if (!should_be_stalled)
                    {
                        for (const Edge &edge : _fwd_edges[rank])
                        {
                            if (_bwdDistances[edge.first] != c::NO_ENTRY && _bwdDistances[edge.first] + edge.second < distance)
                            {
                                should_be_stalled = true;
                                break;
                            }
                        }
                    }
                    if (!should_be_stalled)
                    {
                        for (int edgeid = 0; edgeid < _bwd_shortcuts[rank].size(); edgeid++)
                        {
                            const Shortcut &edge = _bwd_shortcuts[rank][edgeid];
                            Distance new_dist = distance + edge.distance;
                            if (_bwdDistances[edge.node] == c::NO_ENTRY || _bwdDistances[edge.node] > new_dist)
                            {
                                if (_bwdDistances[edge.node] == c::NO_ENTRY)
                                {
                                    bwd_visited_nodes.push_back(edge.node);
                                }
                                _bwdDistances[edge.node] = new_dist;
                                bwd_pQ.push(make_pair(new_dist, edge.node));
                                _added_by_bwd[edge.node] = rank;
                            }
                        }
                        for (int edgeid = 0; edgeid < _bwd_edges[rank].size(); edgeid++)
                        {
                            const Edge &edge = _bwd_edges[rank][edgeid];
                            Distance new_dist = distance + edge.second;
                            if (_bwdDistances[edge.first] == c::NO_ENTRY || _bwdDistances[edge.first] > new_dist)
                            {
                                if (_bwdDistances[edge.first] == c::NO_ENTRY)
                                {
                                    bwd_visited_nodes.push_back(edge.first);
                                }
                                _bwdDistances[edge.first] = new_dist;
                                bwd_pQ.push(make_pair(new_dist, edge.first));
                                _added_by_bwd[edge.first] = rank;
                            }
                        }
                    }
                }
            }
            if (fwd_pQ.empty())
            {
                forward = false;
            }
            else if (bwd_pQ.empty())
            {
                forward = true;
            }
            else
            {
                forward = fwd_pQ.top().first < bwd_pQ.top().first;
            }
        }
        for (auto rank : fwd_visited_nodes)
        {
            _fwdDistances[rank] = c::NO_ENTRY;
        }
        for (auto rank : bwd_visited_nodes)
        {
            _bwdDistances[rank] = c::NO_ENTRY;
        }
        if (midpoint != c::NO_ENTRY)
        {
            struct TmpEdge
            {
                Rank source;
                Rank target;
            };
            shortest_path.reserve(10000);
            vector<TmpEdge> fwd_ch_path;
            fwd_ch_path.reserve(100);
            Rank cur_node = midpoint;
            while (cur_node != source)
            {
                fwd_ch_path.push_back(TmpEdge{
                    _added_by_fwd[cur_node],
                    cur_node,
                });
                cur_node = _added_by_fwd[cur_node];
            }
            vector<TmpEdge> bwd_ch_path;
            bwd_ch_path.reserve(100);
            cur_node = midpoint;
            while (cur_node != target)
            {
                bwd_ch_path.push_back(TmpEdge{
                    cur_node,
                    _added_by_bwd[cur_node],
                });
                cur_node = _added_by_bwd[cur_node];
            }
            vector<TmpEdge> queue;
            queue.reserve(10000);
            for (int i = bwd_ch_path.size() - 1; i >= 0; i--)
            {
                queue.push_back(bwd_ch_path[i]);
            }
            for (int i = 0; i < fwd_ch_path.size(); i++)
            {
                queue.push_back(fwd_ch_path[i]);
            }
            int queue_index = queue.size();
            while (queue_index > 0)
            {
                queue_index--;
                const TmpEdge cur_edge = queue[queue_index];
                bool is_forward = cur_edge.source < cur_edge.target;
                const vector<Shortcut> &shortcuts = is_forward ? _fwd_shortcuts[cur_edge.source] : _bwd_shortcuts[cur_edge.target];
                bool edge_found = false;
                Rank other_side_node = is_forward ? cur_edge.target : cur_edge.source;
                for (const Shortcut &edge : shortcuts)
                {
                    if (edge.node == other_side_node)
                    {
                        edge_found = true;
                        queue[queue_index] = TmpEdge{
                            edge.midnode,
                            cur_edge.target,
                        };
                        queue_index++;
                        if (queue_index == queue.size())
                        {
                            queue.push_back(TmpEdge{
                                cur_edge.source,
                                edge.midnode,
                            });
                        }
                        else
                        {
                            queue[queue_index] = TmpEdge{
                                cur_edge.source,
                                edge.midnode,
                            };
                        }
                        queue_index++;
                        break;
                    }
                }
                if (!edge_found)
                {
                    const vector<Edge> &edges = is_forward ? _fwd_edges[cur_edge.source] : _bwd_edges[cur_edge.target];
                    for (const Edge &edge : edges)
                    {
                        if (edge.first == other_side_node)
                        {
                            edge_found = true;
                            shortest_path.push_back(make_pair(cur_edge.target, edge.second));
                            break;
                        }
                    }
                }
            }
        }
        return shortest_path;
    }

    long nofNodes()
    {
        return _numNodes;
    }

    long nofEdges()
    {
        return _numEdges;
    }

    long nofShortcuts()
    {
        return _numShortcuts;
    }

    long getSizeOfCH()
    {
        return _numEdges * sizeof(Edge) + 4 * sizeof(Rank) * _numNodes + _numShortcuts * sizeof(Shortcut);
    }

    void resetCounter()
    {
        _counter = 0;
    }

    long getCounter()
    {
        return _counter;
    }

    Rank osmIDToRank(long osmID)
    {
        return _osm_id_to_rank[osmID];
    }

private:
    long _numNodes;
    long _numEdges;
    long _numShortcuts;
    long _counter = 0;
    bool _is_initialized = false;
    vector<vector<Edge>> _fwd_edges;
    vector<vector<Edge>> _bwd_edges;
    vector<vector<Shortcut>> _fwd_shortcuts;
    vector<vector<Shortcut>> _bwd_shortcuts;

    vector<Distance> _fwdDistances;
    vector<Distance> _bwdDistances;
    vector<Rank> _added_by_fwd;
    vector<Rank> _added_by_bwd;
    vector<Rank> _osm_id_to_rank;
};

#endif