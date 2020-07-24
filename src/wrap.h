#include<vector>
#include <string>

using namespace std;

// -- Common Geom Vecs
struct geomvecs{
  vector<double> xs;
  vector<double> ys;
};

// we have two basic extraction steps in any solution response. The first is the set
// of nodes. the second is the set of edges. The nodes are typically formed from the
// locations/points attributes. We can extract all these things directly from the pbf
// in a standard format.
// Locations: [id, x, y]

// -- START TSP OBJECTS
struct tspStops{
  void resize(int size){
    sequence.resize(size);
    id.resize(size);
    x.resize(size);
    y.resize(size);
    arrivalTime.resize(size);
    windowStart.resize(size);
    windowEnd.resize(size);
  }
  vector<int> sequence;
  vector<string> id;
  vector<double> x;
  vector<double> y;
  vector<double> arrivalTime;
  vector<double> windowStart;
  vector<double> windowEnd;
};

struct tspEdges{
  void resize(int size){
    geoms.resize(size);
    fromId.resize(size);
    toId.resize(size);
    distance.resize(size);
    time.resize(size);
  }
  vector<geomvecs> geoms;
  vector<string> fromId;
  vector<string> toId;
  vector<double> distance;
  vector<double> time;
};

struct tspTabular{
  tspStops nodes;
  tspEdges edges;
};
// -- END TSP OBJECTS

// -- START IVR OBJECTS
struct stop_block {
  double start = 0;
  double end = 0;
  double slackval = 0;
  double slackcost = 0;
  double tardyval = 0;
  double tardycost = 0;
  double cost = 0;
};
struct edge_block {
  double start = 0;
  double end = 0;
  double cost = 0;
};
struct stop_row {
  int stopid;
  int sequence;
  string location;
  string taskid;
  string jobid;
  string vehicleid;
  int dayIndex = -1;
  string compartmentid;
  unordered_map<string, stop_block> dimVals;
  double x;
  double y;
};
struct edge_row {
  int fromStopId;
  int toStopId;
  geomvecs geom;
  unordered_map<string, edge_block> dimVals;
};
struct ivr_tabular{
  vector<stop_row> route_rows;
  vector<edge_row> edge_rows;
};

struct matrix_tabular{
  vector<string> fids;
  vector<string> tids;
  vector<double> distances;
  vector<double> durations;
};

// -- END IVR OBJECTS

// -- NVD OBJECTS

struct nvd_frontier_item{
  int solutionIndex;
  vector<double> objectiveValues;
  vector<string> objectiveNames;
};
struct nvd_tabular{
  vector<nvd_frontier_item> frontier;
  ivr_tabular tab;
};

// -- END NVD OBJECTS

tspTabular tabulateTSPdata(string& tspSolveRequest, string& solRespString);
tspTabular tabulateTSPTWdata(string& tspSolveRequest, string& solRespString);
ivr_tabular tabulateIVR7(string& ivrSolveRequest, string& solRespString);
ivr_tabular tabulateIVR8(string& ivrSolveRequest, string& solRespString);
nvd_tabular tabulateNVD(string& nvdSolveRequest, string& nvdRespString);
matrix_tabular tabulateMatrix(string& matrixRequest, string& matrixResp);
