#include <Rcpp.h>
#include "wrap.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

double make_precise(double d, double precision) {
  if (precision == 0.0)
    return d;
  if (precision < 0.0) { // round to float, 4-byte precision
    float f = d;
    return (double) f;
  }
  return round(d * precision) / precision;
}

void add_double(std::ostringstream& os, double d, double prec = 0.0) {
  d = make_precise(d, prec); // doubles are ALWAYS coordinates
  const char *cp = (char *)&d;
  os.write((char*) cp, sizeof(double));
}

void add_byte(std::ostringstream& os, char c) {
  os.write((char*) &c, sizeof(char));
}

void add_int(std::ostringstream& os, unsigned int i) {
  const char *cp = (char *)&i;
  os.write((char*) cp, sizeof(int));
}


void write_geomvecs(std::ostringstream& os, geomvecs& data, double prec) {
  auto nrow = data.xs.size();

  add_int(os, nrow);
  for (decltype(nrow) i = 0; i < nrow; i++){
    add_double(os, data.xs[i], prec);
    add_double(os, data.ys[i], prec);
  }
}

void write_data(std::ostringstream& os, geomvecs& sfc, int i = 0, bool EWKB = false,
                int endian = 0, const char *cls = NULL, const char *dim = NULL, double prec = 0.0,
                int srid = 0) {


  add_byte(os, (char) endian);
  unsigned int sf_type = 2;// make_type(cls, dim, EWKB, &tp, srid);
  add_int(os, sf_type);

  write_geomvecs(os, sfc, prec);
}
int native_endian(void) {
  const int one = 1;
  unsigned char *cp = (unsigned char *) &one;
  return (int) *cp;
}

DataFrame geomVecsToSfDf(std::vector<geomvecs>& geoms){
  List output(geoms.size());
  int endian = native_endian();
  bool EWKB = false;
  int srid = 0;
  Function as_sfc("st_as_sfc");
  for (unsigned int i = 0; i < geoms.size(); i++) {
    Rcpp::checkUserInterrupt();
    std::ostringstream os;
    write_data(os, geoms[i], i, EWKB, endian, NULL, NULL, 0.0, srid);
    Rcpp::RawVector raw(os.str().size()); // os -> raw:
    std::string str = os.str();
    const char *cp = str.c_str();
    for (size_t j = 0; j < str.size(); j++){
      raw[j] = cp[j];
    }
    output[i] = raw; // raw vector
  }
  return as_sfc(output);
}


string bytesToString(Vector<RAWSXP>& data){
  char* carray = new char[data.size()];
  for(int i = 0; i < data.size(); i++){
    carray[i] = data[i];
  }
  std::string s(carray, data.size());
  delete[] carray;
  return s;
}

DataFrame tspStopsToDF(tspStops& rows){
  return(DataFrame::create(_["sequence"] = rows.sequence,
                           _["id"] = rows.id,
                           _["x"] = rows.x,
                           _["y"] = rows.y,
                           _["stringsAsFactors"] = false));
}

DataFrame tsptwStopsToDF(tspStops& rows){
  return(DataFrame::create(_["sequence"] = rows.sequence,
                           _["id"] = rows.id,
                           _["x"] = rows.x,
                           _["y"] = rows.y,
                           _["arrivalTimes"] = rows.arrivalTime,
                           _["windowStart"] = rows.windowStart,
                           _["windowEnd"] = rows.windowEnd,
                           _["stringsAsFactors"] = false));
}


DataFrame dropLastRow(DataFrame stopdf){
  // we only have to do this for the tsp stop data frames because the start node
  // and the end node are the same (only in the projection back to the user).
  if(stopdf.nrows() < 1){
    return stopdf;
  }
  IntegerVector rws(stopdf.nrows() - 1);
  for(int i =0 ;i < rws.size(); i++){
    rws[i] = i;
  }
  return (DataFrame::create(
      _["id"] = as<CharacterVector>(stopdf["id"])[rws],
                                                 _["sequence"] = as<NumericVector>(stopdf["sequence"])[rws],
                                                                                                      _["x"] = as<NumericVector>(stopdf["x"])[rws],
                                                                                                                                             _["y"] = as<NumericVector>(stopdf["y"])[rws],
                                                                                                                                                                                    _["stringsAsFactors"] = false
  ));
}

// don't export this function, it's internal
DataFrame tspEdgesToDF(tspStops& stops, tspEdges& edges){
  auto stopdf = tspStopsToDF(stops);
  auto geoms = geomVecsToSfDf(edges.geoms);
  auto edgedf = DataFrame::create(
    _["fromId"] = edges.fromId,
    _["toId"] = edges.toId,
    _["distance"] = edges.distance,
    _["geometry"] = geoms,
    _["stringsAsFactors"] = false);

  Function _left_join("left_join");
  Function _select("select");
  Function _arrange("arrange");
  CharacterVector b1 = CharacterVector::create(_["fromId"] = "id");
  CharacterVector b2 = CharacterVector::create(_["toId"] = "id");
  Function unique("unique");

  auto stopss = dropLastRow(stopdf);

  edgedf = _left_join(edgedf,
                      _select(stopss,
                              _["id"] = "id",
                              _["sequence"] = "sequence",
                              _["fx"] = "x",
                              _["fy"] = "y"), _["by"] = b1);
  edgedf = unique(_left_join(edgedf,
                             _select(stopss,
                                     _["id"] = "id",
                                     _["tx"] = "x",
                                     _["ty"] = "y"), _["by"] = b2));

  return (edgedf);
}


DataFrame tsptwEdgesToDF(tspStops& stops, tspEdges& edges){
  auto stopdf = tspStopsToDF(stops);
  auto geoms = geomVecsToSfDf(edges.geoms);

  auto edgedf = DataFrame::create(
    _["fromId"] = edges.fromId,
    _["toId"] = edges.toId,
    _["distance"] = edges.distance,
    _["time"] = edges.time,
    _["geometry"] = geoms,
    _["stringsAsFactors"] = false);

  Function _left_join("left_join");
  Function _select("select");
  Function _arrange("arrange");
  CharacterVector b1 = CharacterVector::create(_["fromId"] = "id");
  CharacterVector b2 = CharacterVector::create(_["toId"] = "id");

  auto stopss = dropLastRow(stopdf);
  edgedf = _left_join(_["x"] = edgedf,
                      _["y"] = _select(stopss,
                                      _["id"] = "id",
                                      _["sequence"] = "sequence",
                                      _["fx"] = "x",
                                      _["fy"] = "y"),
                                      _["by"] = b1);

  edgedf = _left_join(_["x"] = edgedf,
                      _["y"] = _select(stopss,
                                      _["id"] = "id",
                                      _["tx"] = "x",
                                      _["ty"] = "y"),
                                      _["by"] = b2);
  return (edgedf);
}

DataFrame cleanivrstringfields(DataFrame df, bool compartment){
  Function aschar("as.character");
  df["locationId"] = aschar(df["locationId"]);
  df["taskId"] = aschar(df["taskId"]);
  df["jobId"] = aschar(df["jobId"]);
  df["vehicleId"] = aschar(df["vehicleId"]);
  if(compartment){
    df["compartmentId"] = aschar(df["compartmentId"]);
  }
  return(df);
}

DataFrame ivr7NodesToDF(ivr_tabular& t){
  auto& route_rows = t.route_rows;
  if(route_rows.size() > 0){
    // grab the dimensions we need to iterate over from the first
    // item in the map.
    vector<string> dims;
    for(auto it = route_rows[0].dimVals.begin(); it != route_rows[0].dimVals.end(); it++){
      dims.push_back(it->first);
    }
    // this way we can iterate consistently for all the dimensions.

    // the only bummer is we've collected the data in a row-wise fashion. So we need
    // to convert it into a column-wise   representation.
    NumericVector stopids(route_rows.size());
    NumericVector sequences(route_rows.size());
    StringVector locations(route_rows.size());
    StringVector taskids(route_rows.size());
    StringVector jobids(route_rows.size());
    StringVector vehicleids(route_rows.size());
    NumericVector xs(route_rows.size());
    NumericVector ys(route_rows.size());
    NumericVector dayind(route_rows.size());
    StringVector compartmentids(route_rows.size());
    bool addDayInd = false;
    bool addCompartmentInd = false;

    // index is the dimension, value is
    NumericVector tmp(route_rows.size()); // note to self. In rcpp, this creates a ref vector by default.
    vector<NumericVector> starts(dims.size());
    vector<NumericVector> ends(dims.size());
    vector<NumericVector> slackvals(dims.size());
    vector<NumericVector> slackcosts(dims.size());
    vector<NumericVector> tardyvals(dims.size());
    vector<NumericVector> tardycosts(dims.size());
    vector<NumericVector> costs(dims.size());
    for(unsigned int d = 0; d < dims.size(); d++){
      starts[d] = clone(tmp);
      ends[d] = clone(tmp);
      slackvals[d] = clone(tmp);
      slackcosts[d] = clone(tmp);
      tardyvals[d] = clone(tmp);
      tardycosts[d] = clone(tmp);
      costs[d] = clone(tmp);
    }
    // the primary index is i, the row we're handling. d is the dimension.
    for(unsigned int i = 0; i < route_rows.size(); i++){
      auto& row = route_rows[i];
      stopids(i) = row.stopid;
      sequences(i) = row.sequence;
      locations(i) = row.location;
      taskids(i) = row.taskid;
      jobids(i) = row.jobid;
      vehicleids(i) = row.vehicleid;
      xs(i) = row.x;
      ys(i) = row.y;
      dayind(i) = row.dayIndex;
      if(!addDayInd && row.dayIndex != -1){
        addDayInd = true;
      }
      compartmentids(i) = row.compartmentid;
      if(!addCompartmentInd && row.compartmentid != ""){
        addCompartmentInd = true;
      }

      for(unsigned int d = 0; d < dims.size(); d++){
        auto& sb = row.dimVals[dims[d]];
        starts[d](i) = sb.start;
        ends[d](i) = sb.end;
        slackvals[d](i) = sb.slackval;
        slackcosts[d](i) = sb.slackcost;
        tardyvals[d](i) = sb.tardyval;
        tardycosts[d](i) = sb.tardycost;
        costs[d](i) = sb.cost;
      }
    }

    auto df = DataFrame::create(
      _["stopId"] = stopids,
      _["sequence"] = sequences,
      _["locationId"] = locations,
      _["taskId"] = taskids,
      _["jobId"] = jobids,
      _["vehicleId"] = vehicleids,
      _["x"] = xs,
      _["y"] = ys,
      _["stringsAsFactors"] = false);

    if(addDayInd){
      df.push_back(dayind, "day");
    }
    if(addCompartmentInd){
      df.push_back(compartmentids, "compartmentId");
    }

    for(unsigned int d = 0; d < dims.size(); d++){
      df.push_back(starts[d], (dims[d] + "_start"));
      df.push_back(ends[d], (dims[d] + "_end"));
      df.push_back(slackvals[d], (dims[d] + "_slackval"));
      df.push_back(slackcosts[d], (dims[d] + "_slackcost"));
      df.push_back(tardyvals[d], (dims[d] + "_tardyval"));
      df.push_back(tardycosts[d], (dims[d] + "_tardycost"));
      df.push_back(costs[d], (dims[d] + "_cost"));
    }
    return cleanivrstringfields(df, addCompartmentInd);
  }
  return NULL;
}

DataFrame frontierToDF(vector<nvd_frontier_item>& frontier){
  // we need to translate this into a columnar format.
  vector<NumericVector> objs;
  vector<string> objNames;
  unordered_map<string, int> nLookup;
  NumericVector solnIndex;

  for(auto f : frontier){
    solnIndex.push_back(f.solutionIndex);

    for(unsigned int i = 0; i < f.objectiveNames.size(); i++){
      auto it = nLookup.find(f.objectiveNames[i]);
      int idx = -1;
      if(it == nLookup.end()){
        idx = nLookup.size();
        nLookup.insert(pair<string,int>(f.objectiveNames[i], idx));
        objNames.push_back(f.objectiveNames[i]);
        NumericVector nvec;
        objs.push_back(nvec);
      }else{
        idx = it->second;
      }
      objs[idx].push_back(f.objectiveValues[i]);
    }
  }
  DataFrame df = DataFrame::create(_["SolutionIndex"] = solnIndex);

  for(unsigned int i =0 ; i < objs.size(); i++){
    df[objNames[i]] = objs[i];
  }
  Function asdf("as.data.frame");
  return (asdf(df));
}

DataFrame ivr7EdgesToDF(DataFrame nodes, ivr_tabular& t){
  auto& route_edges = t.edge_rows;
  if(route_edges.size() > 0){
    vector<string> dims;
    for(auto it = route_edges[0].dimVals.begin(); it != route_edges[0].dimVals.end(); it++){
      dims.push_back(it->first);
    }
    vector<geomvecs> allgeoms(route_edges.size());
    for(unsigned int i = 0; i < route_edges.size(); i++){
      allgeoms[i] = route_edges[i].geom;
    }

    NumericVector fromStopIds(route_edges.size());
    NumericVector toStopIds(route_edges.size());

    // index is the dimension, value is
    NumericVector tmp(route_edges.size()); // note to self. In rcpp, this creates a ref vector by default.
    vector<NumericVector> starts(dims.size());
    vector<NumericVector> ends(dims.size());
    vector<NumericVector> costs(dims.size());
    for(unsigned int d = 0; d < dims.size(); d++){
      starts[d] = clone(tmp);
      ends[d] = clone(tmp);
      costs[d] = clone(tmp);
    }
    for(unsigned int i = 0; i < route_edges.size(); i++){
      auto& row = route_edges[i];
      fromStopIds(i) = row.fromStopId;
      toStopIds(i) = row.toStopId;

      for(unsigned int d = 0; d < dims.size(); d++){
        auto& sb = row.dimVals[dims[d]];
        starts[d](i) = sb.start;
        ends[d](i) = sb.end;
        costs[d](i) = sb.cost;
      }
    }
    DataFrame df = DataFrame::create(
      _["fromStopId"] = fromStopIds,
      _["toStopId"] = toStopIds,
      _["geometry"] = geomVecsToSfDf(allgeoms),
      _["stringsAsFactors"] = false);

    for(unsigned int d = 0; d < dims.size(); d++){
      df.push_back(starts[d], (dims[d] + "_start"));
      df.push_back(ends[d], (dims[d] + "_end"));
      df.push_back(costs[d], (dims[d] + "_cost"));
    }

    df = as<DataFrame>(df);

    Function _left_join("left_join");
    Function _select("select");
    CharacterVector b1 = CharacterVector::create(_["fromStopId"] = "stopId");
    CharacterVector b2 = CharacterVector::create(_["toStopId"] = "stopId");
    Function _st_as_sf("st_as_sf");
    df = _st_as_sf(df); // cooerce this out of a "list" format.

    auto colnames = as<vector<string>>(nodes.names());
    bool hasDayInd = false;
    for(auto cn : colnames){
      if(cn == "day"){
        hasDayInd = true;
        break;
      }
    }
    Rcout << endl;
    if(hasDayInd){
      df = _left_join(_["x"] = df,
                      _["y"] = _select(nodes,
                                      _["stopId"] = "stopId",
                                      _["sequence"] = "sequence",
                                      _["vehicleId"] = "vehicleId",
                                      _["fx"] = "x",
                                      _["fy"] = "y",
                                      _["fromLocationId"] = "locationId",
                                      _["day"] = "day"),
                                      _["by"] = b1);
    }else{
      df = _left_join(_["x"] = df,
                      _["y"] = _select(nodes,
                                      _["stopId"] = "stopId",
                                      _["sequence"] = "sequence",
                                      _["vehicleId"] = "vehicleId",
                                      _["fx"] = "x",
                                      _["fy"] = "y",
                                      _["fromLocationId"] = "locationId"),
                                      _["by"] = b1);
    }

    df = _left_join(_["x"] = df,
                    _["y"] = _select(nodes,
                                    _["stopId"] = "stopId",
                                    _["tx"] = "x",
                                    _["ty"] = "y",
                                    _["toLocationId"] = "locationId"),
                                    _["by"] = b2);

    return(df);
  }
  return (NULL);
}

DataFrame matrixToDf(string& indata, string& outdata){
  matrix_tabular tab = tabulateMatrix(indata, outdata);
  return DataFrame::create( Named("fromId") = tab.fids,
                            Named("toId") = tab.tids,
                            Named("distance") = tab.distances,
                            Named("duration") = tab.durations,
                            _["stringsAsFactors"] = false);
}


//' @name tabulate
//' @title tabulates a solution response (requires the base request)
//' @param solResponse The solution response from the api.
//' @param solveRequest The base-model which was solved (i.e. submitted to the api).
//' @export
// [[Rcpp::export]]
List tabulate(Reference solResponse, Reference solveRequest){
  std::string srType = as<std::string>(solveRequest.slot("type"));
  Function tobytes("toBytes");
  Vector<RAWSXP> solResp = tobytes(solResponse);
  Vector<RAWSXP> sr = tobytes(solveRequest);
  Function _tab_transitrules("transitRulesToTable");
  Function _tab_infeasibilities("infeasibilitiesToTable");
  if(srType == "TSP.SolveRequest"){
    string outdata = bytesToString(solResp);
    string indata =  bytesToString(sr);

    auto t = tabulateTSPdata(indata, outdata);
    return(List::create( _["edges"] = tspEdgesToDF(t.nodes, t.edges),
                         _["nodes"] = tspStopsToDF(t.nodes)));
  }
  if(srType == "TSPTW.SolveRequest"){
    string outdata = bytesToString(solResp);
    string indata =  bytesToString(sr);
    auto t = tabulateTSPTWdata(indata, outdata);
    return(List::create( _["edges"] = tsptwEdgesToDF(t.nodes, t.edges),
                         _["nodes"] = tsptwStopsToDF(t.nodes)));
  }
  if(srType == "IVR7.SolveRequest"){
    string outdata = bytesToString(solResp);
    string indata =  bytesToString(sr);
    auto t = tabulateIVR7(indata, outdata);
    auto ns = ivr7NodesToDF(t);
    return List::create(_["nodes"] = ns,
                        _["edges"] = ivr7EdgesToDF(ns, t),
                        _["transitRules"] = _tab_transitrules(solResponse),
                        _["infeasibilities"] = _tab_infeasibilities(solResponse));
  }
  if(srType == "IVR8.SolveRequest"){
    string outdata = bytesToString(solResp);
    string indata =  bytesToString(sr);
    auto t = tabulateIVR8(indata, outdata);
    auto ns = ivr7NodesToDF(t); // this is okay, it knows how to do both.
    Function _tab_compartments("compartmentsToTable");
    return List::create(_["nodes"] = ns,
                        _["edges"] = ivr7EdgesToDF(ns, t),
                        _["transitRules"] = _tab_transitrules(solResponse),
                        _["compartmentSummary"] = _tab_compartments(ns, solResponse, solveRequest),
                        _["infeasibilities"] = _tab_infeasibilities(solResponse));
  }
  if(srType == "CVRP.SolveRequest"){
    Function tabcvrp("tabulatecvrp");
    return(tabcvrp(_["sr"] = solveRequest, _["resp"] = solResponse));
  }
  if(srType == "NVD.SolveRequest"){
    string outdata = bytesToString(solResp);
    string indata =  bytesToString(sr);
    auto t = tabulateNVD(indata, outdata);
    if(t.frontier.size() == 0){
      auto ns = ivr7NodesToDF(t.tab);
      return List::create(_["nodes"] = ns,
                          _["edges"] = ivr7EdgesToDF(ns, t.tab));
    }else{
      return(List::create(_["frontier"] = frontierToDF(t.frontier)));
    }
  }
  if(srType == "Matrix.MatrixRequest"){
    string outdata = bytesToString(solResp);
    string indata =  bytesToString(sr);
    auto t = matrixToDf(indata, outdata); // long form tabulation
    return List::create(_["elements"] = t);
  }
  return NULL;
}


