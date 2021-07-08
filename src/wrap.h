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

// Frontier objects
struct frontier_item{
  int solutionIndex;
  vector<double> objectiveValues;
  vector<string> objectiveNames;
};

// -- NVD OBJECTS
struct nvd_tabular{
  vector<frontier_item> frontier;
  ivr_tabular tab;
};
// -- END NVD OBJECTS


// -- NDD OBJECTS
struct ndd_tabular{
  vector<frontier_item> frontier;
  ivr_tabular tab;
};
// -- END NDD OBJECTS

// -- NS3 OBJECTS
// just a note, R likes things in columnar format.
// this this effectively extracts everything in columnar format
// so that we can just coerce it into a data-frame afterwards
struct ns3_nodeflow_tab{
  vector<string> nodeId;
  vector<double> inFlow;
  vector<double> outFlow;
  vector<double> flowCost;
  vector<double> fixedCost;
  vector<double> productFlowCost;
  vector<double> productFixedCost;
  vector<double> productionAmount;
  vector<double> productionPenalty;
  vector<double> productionCost;
  vector<double> consumptionAmount;
  vector<double> consumptionPenalty;
  vector<double> consumptionCost;
};
struct ns3_nodepflow_tab{
  vector<string> nodeId;
  vector<string> productId;
  vector<float> inFlow;
  vector<float> outFlow;
  vector<float> flowCost;
  vector<float> fixedCost;
  vector<float> productionAmount;
  vector<float> productionPenalty;
  vector<float> productionCost;
  vector<float> consumptionAmount;
  vector<float> consumptionPenalty;
  vector<float> consumptionCost;
};
struct ns3_assignmnet_tab{
  vector<string> source;
  vector<string> destination;
  vector<string> productId;
  vector<double> amount;
  vector<double> cost;
  vector<string> laneRateId;
  vector<string> costModelId;
  vector<double> distance;
  vector<double> duration;
};
struct ns3_geom_tab{
  vector<string> fromId;
  vector<string> toId;
  vector<geomvecs> geom;
};
struct ns3_prod_transform_assigments{
  vector<string> nodeId;
  vector<string> productTransformId;
  vector<string> productId;
  vector<string> transformState; //either input or output.
  vector<double> amount;
  vector<double> cost;
  vector<double> fixedCost;
  vector<double> penaltyAmount;
  vector<double> penaltyCost;
};
struct ns3_nodes{
  vector<string> Id;
  vector<double> xs;
  vector<double> ys;
};
struct ns3_tabular{
  ns3_nodes nodes;
  ns3_nodeflow_tab nodeflowtab;
  ns3_nodepflow_tab nodeflowproducttab;
  ns3_assignmnet_tab assignmenttab;
  ns3_geom_tab geomtab;
  ns3_prod_transform_assigments prodtransformtab;
};
// -- END NS3 OBJECTS

tspTabular tabulateTSPdata(string& tspSolveRequest, string& solRespString);
tspTabular tabulateTSPTWdata(string& tspSolveRequest, string& solRespString);
ivr_tabular tabulateIVR7(string& ivrSolveRequest, string& solRespString);
ivr_tabular tabulateIVR8(string& ivrSolveRequest, string& solRespString);
nvd_tabular tabulateNVD(string& nvdSolveRequest, string& nvdRespString);
ndd_tabular tabulateNDD(string& nddSolveRequest, string& nddRespString);
matrix_tabular tabulateMatrix(string& matrixRequest, string& matrixResp);
ns3_tabular tabulateNS3(string& ns3Request, string& ns3Resp);
