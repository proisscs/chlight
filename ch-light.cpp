#include <iostream>
#include "ch-light.h"
#include "dijkstra.h"
#include "ch.h"
#include "CHParser.h"
#include <chrono>
#include <ctime>
#include <random>
#include <sstream>

using namespace std;

string getTime()
{
    auto time1 = std::chrono::system_clock::now();
    time_t time2 = chrono::system_clock::to_time_t(time1);
    return ctime(&time2);
}

string getQueryTypeString(int query_type)
{
    if (query_type == c::QUERY_CHLIGHTDIST)
    {
        return "CHLight distance only";
    }
    else if (query_type == c::QUERY_CHLIGHT)
    {
        return "CHLight";
    }
    else if (query_type == c::QUERY_DIJK)
    {
        return "Bidirectional Dijkstra";
    }
    else if (query_type == c::QUERY_CHLIGHTSOD)
    {
        return "LCH with stall-on-demand";
    }
    else if (query_type == c::QUERY_UNIDIJK)
    {
        return "Normal Dijkstra";
    }
    else if (query_type == c::QUERY_UNIDIJKDIST)
    {
        return "Normal Dijkstra Distance only";
    }
    else if (query_type == c::QUERY_BACKUNIDIJK)
    {
        return "Backward Dijkstra";
    }
    else if (query_type == c::QUERY_CH)
    {
        return "CH";
    }
    else if (query_type == c::QUERY_CHDIST)
    {
        return "CH Distance only";
    }
    return "Unknown";
}

class Tester
{
public:
    void constructCHLight(string input)
    {
        _chlight.constructCHLight(input);
        _num_nodes = _chlight.nofNodes();
    }
    void readInDijkstra(string input)
    {
        _dijk.readFromFMIFile(input);
        _num_nodes = _dijk.nofNodes();
    }
    void readInCH(string input)
    {
        _ch.constructCH(input);
        _num_nodes = _ch.nofNodes();
    }
    void generateRandomQueries(int num_queries)
    {
        if (_num_nodes == c::NO_ENTRY)
        {
            cout << "No data loaded yet, abort..." << endl;
            return;
        }
        srand(time(nullptr));
        _sources.resize(num_queries);
        _targets.resize(num_queries);
        for (int i = 0; i < num_queries; i++)
        {
            int source = rand() % _num_nodes;
            int target = rand() % _num_nodes;
            _sources[i] = source;
            _targets[i] = target;
        }
    }
    void readQueriesFromFile(string input)
    {
        if (_num_nodes == c::NO_ENTRY)
        {
            cout << "No data loaded yet, abort..." << endl;
            return;
        }
        srand(time(nullptr));
        _sources.clear();
        _targets.clear();
        Rank source, target;
        ifstream reader(input);
        while (reader >> source >> target)
        {
            _sources.push_back(source);
            _targets.push_back(target);
        }
        reader.close();
    }

    void writeQueriesToFile(string output)
    {
        ofstream writer(output);
        for (int i = 0; i < _sources.size(); i++)
        {
            writer << _sources.at(i) << " " << _targets.at(i) << "\n";
        }
        writer.close();
    }

    pair<double, long> speedTest(int query_type = c::QUERY_CHLIGHTDIST)
    {
        if (!_checkInitialization(query_type))
        {
            return make_pair(-1, -1);
        }
        Timer timer;
        vector<int> results(_targets.size());
        vector<pair<Rank, Rank>> queries = _getQueriesForType(query_type);
        if (query_type == c::QUERY_CHLIGHTDIST)
        {
            timer.start();
            for (int i = 0; i < _sources.size(); i++)
            {
                results[i] = _chlight.getDistance(queries[i].first, queries[i].second);
            }
            timer.stop();
        }
        else if (query_type == c::QUERY_CHLIGHT)
        {
            timer.start();
            for (int i = 0; i < _sources.size(); i++)
            {
                results[i] = _chlight.getShortestPath(queries[i].first, queries[i].second).size();
            }
            timer.stop();
        }
        else if (query_type == c::QUERY_CHLIGHTSOD)
        {
            timer.start();
            for (int i = 0; i < _sources.size(); i++)
            {
                results[i] = _chlight.getShortestPathWithSoD(queries[i].first, queries[i].second).size();
            }
            timer.stop();
        }
        else if (query_type == c::QUERY_DIJK)
        {
            timer.start();
            for (int i = 0; i < _sources.size(); i++)
            {
                results[i] = _dijk.getShortestPathDijkstra(queries[i].first, queries[i].second).size();
            }
            timer.stop();
        }
        else if (query_type == c::QUERY_UNIDIJK)
        {
            timer.start();
            for (int i = 0; i < _sources.size(); i++)
            {
                results[i] = _dijk.getShortestPathUnidirectionalDijkstra(queries[i].first, queries[i].second).size();
            }
            timer.stop();
        }
        else if (query_type == c::QUERY_UNIDIJKDIST)
        {
            timer.start();
            for (int i = 0; i < _sources.size(); i++)
            {
                results[i] = _dijk.getDistanceUnidirectionalDijkstra(queries[i].first, queries[i].second);
            }
            timer.stop();
        }
        else if (query_type == c::QUERY_CHDIST)
        {
            timer.start();
            for (int i = 0; i < _sources.size(); i++)
            {
                results[i] = _ch.getDistance(queries[i].first, queries[i].second);
            }
            timer.stop();
        }
        else if (query_type == c::QUERY_CH)
        {
            timer.start();
            for (int i = 0; i < _sources.size(); i++)
            {
                results[i] = _ch.getShortestPath(queries[i].first, queries[i].second).size();
            }
            timer.stop();
        }
        else if (query_type == c::QUERY_DIJKALL)
        {
            timer.start();
            for (int i = 0; i < _sources.size(); i++)
            {
                results[i] = _dijk.oneToAllDijkstra(queries[i].first);
            }
            timer.stop();
        }
        long sumDistances = 0;
        for (int i = 0; i < _sources.size(); i++)
        {
            sumDistances += results[i];
        }
        return make_pair(timer.secs() / _sources.size(), sumDistances);
    }

    vector<long> getVisitedNodes(int query_type)
    {
        vector<long> result;
        if (!_checkInitialization(query_type))
        {
            return result;
        }
        result.resize(_sources.size());

        for (int i = 0; i < _sources.size(); i++)
        {
            long res = -1;
            if (query_type == c::QUERY_CH)
            {
                res = _ch.getNumNodesInUnpackedSearchTreeAndExpandedNodes(_sources[i], _targets[i]).second;
            }
            else if (query_type == c::QUERY_DIJK)
            {
                res = _dijk.getNumberOfVisitedNodes(_sources[i], _targets[i]);
            }
            else if (query_type == c::QUERY_CHLIGHT)
            {
                res = _chlight.getNumberVisitedNodes(_sources[i], _targets[i]);
            }
            else if (query_type == c::QUERY_CHDIST)
            {
                res = _ch.getNumberVisitedNodes(_sources[i], _targets[i]);
            }
            result[i] = res;
        }
        return result;
    }

    long getSize(int query_type)
    {
        if (!_checkInitialization(query_type))
        {
            return -1;
        }
        if (query_type == c::QUERY_CHLIGHTDIST)
        {
            return _chlight.getSizeOfCH();
        }
        else if (query_type == c::QUERY_CHLIGHT)
        {
            return _chlight.getSizeOfCH();
        }
        else if (query_type == c::QUERY_CHLIGHTSOD)
        {
            return _chlight.getSizeOfCH();
        }
        else if (query_type == c::QUERY_DIJK)
        {
            return _dijk.getSizeOfDijkstra();
        }
        else if (query_type == c::QUERY_UNIDIJK)
        {
            return _dijk.getSizeOfUniDijkstra();
        }
        else if (query_type == c::QUERY_UNIDIJKDIST)
        {
            return _dijk.getSizeOfUniDijkstra();
        }
        else if (query_type == c::QUERY_CH)
        {
            return _ch.getSizeOfCH();
        }
        else if (query_type == c::QUERY_CHDIST)
        {
            return _ch.getSizeOfCH();
        }
        return -1;
    }
    long getCounter(int query_type)
    {
        if (!_checkInitialization(query_type))
        {
            return -1;
        }
        if (query_type == c::QUERY_CHLIGHTDIST || query_type == c::QUERY_CHLIGHT)
        {
            return _chlight.getCounter();
        }
        if (query_type == c::QUERY_CH || query_type == c::QUERY_CHDIST)
        {
            return _ch.getCounter();
        }
        return _dijk.getCounter();
    }
    void resetCounter(int query_type)
    {
        if (!_checkInitialization(query_type))
        {
            return;
        }
        if (query_type == c::QUERY_CHLIGHTDIST || query_type == c::QUERY_CHLIGHT)
        {
            _chlight.resetCounter();
        }
        else if (query_type == c::QUERY_CH || query_type == c::QUERY_CHDIST)
        {
            _ch.resetCounter();
        }
        else
        {
            _dijk.resetCounter();
        }
    }
    long getNumEdges(int query_type)
    {
        if (!_checkInitialization(query_type))
        {
            return -1;
        }
        if (query_type == c::QUERY_CHLIGHTDIST || query_type == c::QUERY_CHLIGHT)
        {
            return _chlight.nofEdges();
        }
        if (query_type == c::QUERY_CH || query_type == c::QUERY_CHDIST)
        {
            return _ch.nofEdges() + _ch.nofShortcuts();
        }
        if (query_type == c::QUERY_DIJK)
        {
            return 2 * _dijk.nofEdges();
        }
        return _dijk.nofEdges();
    }
    string test(int query_type)
    {
        if (!_checkInitialization(query_type))
        {
            return "";
        }
        resetCounter(query_type);
        double average_query_time = speedTest(query_type).first;
        double size = getSize(query_type);
        string output = "";
        output += "Query type: " + getQueryTypeString(query_type) + "\n";
        output += "Average query time: " + to_string(average_query_time) + " seconds\n";
        output += "Num edges: " + to_string(getNumEdges(query_type)) + "\n";
        output += "Size: " + to_string(size / 1024. / 1024. / 1024.) + " GB\n";
        output += "Average number settled nodes: " + to_string(getCounter(query_type) / _sources.size()) + "\n";
        return output;
    }
    bool debugWithFMI(int query_type)
    {
        vector<int> results(_sources.size());
        vector<int> results_fmi(_sources.size());
        if (!_ch_fmi.isInitialized())
        {
            cout << "CHFMI not initialized, abort..." << endl;
            return false;
        }
        if (!_checkInitialization(query_type))
        {
            return false;
        }
        vector<pair<Rank, Rank>> queries = _getQueriesForType(query_type);
        for (int i = 0; i < _sources.size(); i++)
        {
            if (query_type == c::QUERY_CHLIGHTDIST)
            {
                results[i] = _chlight.getDistance(queries[i].first, queries[i].second);
            }
            else if (query_type == c::QUERY_CHLIGHT)
            {
                vector<Edge> shortest_path = _chlight.getShortestPath(queries[i].first, queries[i].second);
                if (shortest_path.size() == 0)
                {
                    results[i] = c::NO_ENTRY;
                }
                else
                {
                    Distance dist = 0;
                    for (int j = 0; j < shortest_path.size(); j++)
                    {
                        dist += shortest_path[j].second;
                    }
                    results[i] = dist;
                }
            }
            else if (query_type == c::QUERY_CHLIGHTSOD)
            {
                vector<Edge> shortest_path = _chlight.getShortestPathWithSoD(queries[i].first, queries[i].second);
                if (shortest_path.size() == 0)
                {
                    results[i] = c::NO_ENTRY;
                }
                else
                {
                    Distance dist = 0;
                    for (int j = 0; j < shortest_path.size(); j++)
                    {
                        dist += shortest_path[j].second;
                    }
                    results[i] = dist;
                }
            }
            else if (query_type == c::QUERY_DIJK)
            {
                vector<Edge> shortest_path = _dijk.getShortestPathDijkstra(queries[i].first, queries[i].second);
                if (shortest_path.size() == 0)
                {
                    results[i] = c::NO_ENTRY;
                }
                else
                {
                    Distance dist = 0;
                    for (int j = 0; j < shortest_path.size(); j++)
                    {
                        dist += shortest_path[j].second;
                    }
                    results[i] = dist;
                }
            }
            else if (query_type == c::QUERY_UNIDIJK)
            {
                vector<Edge> shortest_path = _dijk.getShortestPathUnidirectionalDijkstra(queries[i].first, queries[i].second);
                if (shortest_path.size() == 0)
                {
                    results[i] = c::NO_ENTRY;
                }
                else
                {
                    Distance dist = 0;
                    for (int j = 0; j < shortest_path.size(); j++)
                    {
                        dist += shortest_path[j].second;
                    }
                    results[i] = dist;
                }
            }
            else if (query_type == c::QUERY_CH)
            {
                vector<Edge> shortest_path = _ch.getShortestPath(queries[i].first, queries[i].second);
                if (shortest_path.size() == 0)
                {
                    results[i] = c::NO_ENTRY;
                }
                else
                {
                    Distance dist = 0;
                    for (int j = 0; j < shortest_path.size(); j++)
                    {
                        dist += shortest_path[j].second;
                    }
                    results[i] = dist;
                }
            }
            else if (query_type == c::QUERY_UNIDIJKDIST)
            {
                results[i] = _dijk.getDistanceUnidirectionalDijkstra(queries[i].first, queries[i].second);
            }
            else if (query_type == c::QUERY_CHDIST)
            {
                results[i] = _ch.getDistance(_sources[i], _targets[i]);
            }
            results_fmi[i] = _ch_fmi.getDistance(_ch_fmi.osmIDToRank(_sources[i]), _ch_fmi.osmIDToRank(_targets[i]));
        }
        for (int i = 0; i < _sources.size(); i++)
        {
            if (results[i] != results_fmi[i])
            {
                cout << "Error at query " << i + 1 << ": " << _sources[i] << " " << _targets[i] << " " << results[i] << " " << results_fmi[i] << " " << results[i] - results_fmi[i] << endl;
                return false;
            }
            else
            {
                cout << "Correct: " << queries[i].first << " " << queries[i].second << " " << results[i] << endl;
            }
        }
        return true;
    }

    void createDijkstraRankQueryFile(string output_prefix, int num_queries_per_rank = 100)
    {
        _checkInitialization(c::QUERY_UNIDIJKDIST, false);
        long rank = 64;
        vector<long> dijkstra_ranks;
        dijkstra_ranks.push_back(rank);
        long limit = 3 * _dijk.nofNodes() / 7;
        while (rank < limit)
        {
            rank *= 2;
            dijkstra_ranks.push_back(rank);
        }
        vector<vector<pair<Rank, Rank>>> queries_per_rank(dijkstra_ranks.size());
        for (int i = 0; i < dijkstra_ranks.size(); i++)
        {
            queries_per_rank[i].reserve(num_queries_per_rank);
        }
        for (int i = 0; i < num_queries_per_rank;)
        {
            Rank source = rand() % _dijk.nofNodes();
            vector<Rank> settled_nodes = _dijk.getNodesOfDijkstraRanks(source, dijkstra_ranks);
            if (settled_nodes[settled_nodes.size() - 1] != c::NO_ENTRY)
            {
                i++;
                for (int j = 0; j < settled_nodes.size(); j++)
                {
                    queries_per_rank[j].push_back(make_pair(source, settled_nodes[j]));
                }
            }
        }
        for (int i = 0; i < queries_per_rank.size(); i++)
        {
            ofstream writer(output_prefix + "-rank" + to_string(i + 6) + ".txt");
            for (int j = 0; j < queries_per_rank[i].size(); j++)
            {
                writer << queries_per_rank[i][j].first << " " << queries_per_rank[i][j].second << "\n";
            }
            writer.close();
        }
    }

    void printGraphStats(int query_type = c::QUERY_CH)
    {
        if (!_checkInitialization(query_type, false))
        {
            return;
        }
        if (query_type == c::QUERY_CH)
        {
            cout << "Number of nodes: " << _ch.nofNodes() << endl;
            cout << "Number of edges: " << _ch.nofEdges() << endl;
            cout << "Number of shortcuts: " << _ch.nofShortcuts() << endl;
        }
        else if (query_type == c::QUERY_DIJK)
        {
            cout << "Number of nodes: " << _dijk.nofNodes() << endl;
            cout << "Number of edges: " << _dijk.nofEdges() << endl;
        }
        else if (query_type == c::QUERY_CHLIGHT)
        {
            cout << "Number of nodes: " << _chlight.nofNodes() << endl;
            cout << "Number of edges: " << _chlight.nofEdges() << endl;
        }
    }

    void printGraphSize()
    {
        cout << "Sizes in Bytes/GB:" << endl;
        long divisor = 1024 * 1024 * 1024;
        if (_checkInitialization(c::QUERY_CH, false, false))
        {
            cout << "CH: " << _ch.getSizeOfCH() << " " << (double)(_ch.getSizeOfCH()) / divisor << endl;
        }
        if (_checkInitialization(c::QUERY_DIJK, false, false))
        {
            cout << "BiDijk: " << _dijk.getSizeOfDijkstra() << " " << (double)(_dijk.getSizeOfDijkstra()) / divisor << endl;
            cout << "UniDijk: " << _dijk.getSizeOfUniDijkstra() << " " << (double)(_dijk.getSizeOfUniDijkstra()) / divisor << endl;
        }
        if (_checkInitialization(c::QUERY_CHLIGHT, false, false))
        {
            cout << "LCH: " << _chlight.getSizeOfCH() << " " << (double)(_chlight.getSizeOfCH()) / divisor << endl;
        }
    }

    void writeUpwardSearchToFile(string file_name, Rank source)
    {
        _checkInitialization(c::QUERY_CHLIGHT, false, false);
        vector<pair<Rank, Rank>> upward_search = _chlight.getUpwardSearch(source);
        ofstream writer(file_name);
        for (const auto &pair : upward_search)
        {
            writer << pair.first << " " << pair.second << "\n";
        }
        writer.close();
    }

    void getMaxLevel()
    {
        cout << "Max level: " << _ch_fmi.maxLevel << endl;
        if(_ch_fmi.maxLevel > 255)
        {
            long count_nodes = 0;
            for(long i=0; i<_ch_fmi.nofNodes(); i++)
            {
                if(_ch_fmi.nodeIDToLevel[i] > 255)
                {
                    count_nodes++;
                }
            }
            cout << "Nodes above level 255: " << count_nodes << endl;
        }
    }

private:
    CHLight _chlight;
    Dijkstra _dijk;
    CH _ch;
    vector<int> _sources;
    vector<int> _targets;
    CHParser _ch_fmi;
    long _num_nodes = c::NO_ENTRY;

    bool _checkInitialization(int query_type, bool check_queries = true, bool verbose = true)
    {
        if (_sources.size() == 0 && check_queries)
        {
            if (verbose)
            {
                cout << "No queries available, abort..." << endl;
            }
            return false;
        }
        if ((query_type == c::QUERY_CHLIGHTDIST || query_type == c::QUERY_CHLIGHT || query_type == c::QUERY_CHLIGHTSOD) && !_chlight.isInitialized())
        {
            if (verbose)
            {
                cout << "CHLight not initialized, abort..." << endl;
            }
            return false;
        }
        if ((query_type == c::QUERY_CH || query_type == c::QUERY_CHDIST) && !_ch.isInitialized())
        {
            if (verbose)
            {
                cout << "CH not initialized, abort..." << endl;
            }
            return false;
        }
        if ((query_type == c::QUERY_DIJK || query_type == c::QUERY_UNIDIJK || query_type == c::QUERY_UNIDIJKDIST) && !_dijk.isInitialized())
        {
            if (verbose)
            {
                cout << "Dijkstra not initialized, abort..." << endl;
            }
            return false;
        }
        return true;
    }

    vector<pair<Rank, Rank>> _getQueriesForType(int query_type)
    {
        vector<pair<Rank, Rank>> queries(_sources.size());
        for (int i = 0; i < _sources.size(); i++)
        {
            Rank source, target;
            if (query_type == c::QUERY_CH)
            {
                source = _ch.osmIDToRank(_sources[i]);
                target = _ch.osmIDToRank(_targets[i]);
            }
            else if (query_type == c::QUERY_DIJK || query_type == c::QUERY_DIJKALL || query_type == c::QUERY_UNIDIJK || query_type == c::QUERY_UNIDIJKDIST)
            {
                source = _dijk.osmIDToRank(_sources[i]);
                target = _dijk.osmIDToRank(_targets[i]);
            }
            else if (query_type == c::QUERY_CHLIGHT || query_type == c::QUERY_CHLIGHTDIST || query_type == c::QUERY_CHLIGHTSOD)
            {
                source = _chlight.osmIDToRank(_sources[i]);
                target = _chlight.osmIDToRank(_targets[i]);
            }
            queries[i] = make_pair(source, target);
        }
        return queries;
    }
};

int main(int argc, char *argv[])
{
    string input;
    string input_chfmi = "none", input_dijk = "none";
    int mode = 0;
    int num_runs = 100;
    int kept_shortcuts = 0;
    uint seed = time(nullptr);

    cout << "Usage: " << argv[0] << " <graph> <ch>" << endl;
    input_chfmi = argv[1];
    input_dijk = argv[2];
    Tester tester;
    tester.constructCHLight(input_chfmi);
    tester.readInCH(input_chfmi);
    tester.readInDijkstra(input_dijk);
    tester.generateRandomQueries(num_runs);
    pair<double, long> result_ch = tester.speedTest(c::QUERY_CH);
    pair<double, long> result_chl = tester.speedTest(c::QUERY_CHLIGHT);
    pair<double, long> result_dijk = tester.speedTest(c::QUERY_DIJK);
    cout << "Run time CH: " << result_ch.first << ", checksum " << result_ch.second << endl;
    cout << "Run time CHLight: " << result_chl.first << ", checksum " << result_chl.second << endl;
    cout << "Run time BiDijk: " << result_dijk.first << ", checksum " << result_dijk.second << endl;
    tester.printGraphSize();
}
