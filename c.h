#ifndef FILE_C_SEEN
#define FILE_C_SEEN

#include <climits>
#include <boost/archive/binary_iarchive.hpp>
#include <vector>

namespace c
{
    int NO_ENTRY = -1;
    int QUERY_CHLIGHTDIST = 0;
    int QUERY_CHLIGHT = 1;
    int QUERY_DIJK = 2;
    int QUERY_CHLIGHTSOD = 3;
    int QUERY_UNIDIJK = 4;
    int QUERY_UNIDIJKDIST = 5;
    int QUERY_BACKUNIDIJK = 6;
    int QUERY_CH = 7;
    int QUERY_CHDIST = 8;
    int QUERY_DIJKALL = 9;
}

typedef int Rank;
typedef int NodeID;
typedef int EdgeID;
typedef int Distance;
typedef unsigned char Level;
typedef int Group;

struct Shortcut
{
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar &node;
        ar &distance;
        ar &midnode;
    }
    Rank node;
    Distance distance;
    Rank midnode;
};

typedef std::pair<Rank, Distance> Edge;

std::vector<Rank> getOSMIDToRankMapping(std::vector<std::pair<long, Rank>> &osm_rank_pairs)
{
    long num_nodes = osm_rank_pairs.size();
    sort(osm_rank_pairs.begin(), osm_rank_pairs.end());
    std::vector<Rank> osm_id_to_rank(num_nodes);
    for (long i = 0; i < num_nodes; i++)
    {
        osm_id_to_rank[i] = osm_rank_pairs[i].second;
    }
    return osm_id_to_rank;
}

#endif