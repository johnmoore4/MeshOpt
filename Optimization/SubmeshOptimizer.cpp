#include "SubmeshOptimizer.h"
#include "SubmeshMeritEvaluator.h"
#include "LocalMeshMeritEvaluator.h"
#include "Mel.h"
#include "MeshContainer.h"
#include "LBFGSOptimizer.h"

using namespace std;

int SubmeshOptimizer::RandomizeNodes(double factor){
  
  std::unordered_set<const MEl*> watch;
  std::map<gind,gind> activenodes;
  
 
  node_map& nodes = mesh.getNodesNC();
  //const ideal_map& idealElements = mesh.getIdealElements();
  //std::map<MEl*,double>& all_merits = mesh.getAllMeritsNC();

  FindActiveElements(watch,activenodes,1.3,2,-1);
  for(auto nd = activenodes.begin(); nd != activenodes.end(); ++nd){
    arma::vec3 xyz = nodes.at(nd->first)->getXYZ();
    for(int i=0; i<2; i++){
      xyz(i)+= factor*((double) std::rand() / (RAND_MAX));
    }
    nodes.at(nd->first)->setXYZ(xyz.memptr());
  }

}

int SubmeshOptimizer::OptimizeGlobal(double threshold, int Nnearest, 
				     int type, double factor){

  std::unordered_set<const MEl*> watch;
  std::map<gind,gind> activenodes;
  //const ideal_map& idealElements = mesh.getIdealElements();
  //std::map<MEl*,double>& all_merits = mesh.getAllMeritsNC();

  FindActiveElements(watch,activenodes,threshold,Nnearest,type);

  cout << "watch size: " << watch.size() << endl;
  cout << "active size: " << activenodes.size() << endl;

  std::vector<const MEl*> temp(watch.size());
  int cnt = 0;
  for(auto el = watch.begin(); el != watch.end(); ++el, ++cnt){
    temp[cnt] = *el;
  }
  LocalMeshMeritEvaluator lm_evaluator(activenodes,temp,
				       idealElements,optel_manager,factor);
  
  LBFGSOptimizer global_optimizer(lm_evaluator);
  global_optimizer.plot_exit_flag=true;
  global_optimizer.Optimize();
  {
    int cnt = 0;
    const std::vector<double>& merits = lm_evaluator.getMerits();
    for(auto el = watch.begin(); el != watch.end(); ++el, ++cnt){
      all_merits[*el] = merits[cnt];
    }
  }
}

double SubmeshOptimizer::OptimizeOneByOne(double threshold, int Nnearest,
				       int type, double factor){

  std::unordered_set<const MEl*> watch;
  std::map<gind,gind> activenodes;
  
 
  const node_map& nodes = mesh.getNodes();
  //const ideal_map& idealElements = mesh.getIdealElements();
  //std::map<MEl*,double>& all_merits = mesh.getAllMeritsNC();
  

  arma::wall_clock timer;
  timer.tic();
  FindActiveElements(watch,activenodes,threshold,Nnearest,type);


  std::map<gind, std::vector<const MEl*> > Node2ElementList = 
    CreateNode2ElementList(watch,activenodes);
 
  std::cout << "active nodes size: " << activenodes.size() << std::endl;

  for(auto nd = activenodes.begin(); nd != activenodes.end(); ++nd){

    const std::vector<const MEl*>& elements = Node2ElementList[nd->first];
    const int sz = elements.size();

    if(sz == 1 && nodes.at(nd->first)->getType() == 2){
      //std::cout << "size is 1!" << std::endl;
      //std::cout << nodes.at(nd->first)->getGeoEntity() << std::endl;
    }
    bool to_opt = true;
    //if(elements[0]->getOrder() > 1 && sz != 1) to_opt = false;

    if(to_opt){
    //if(sz == 2){
    std::vector<ElIdealPair> temp;
    for(int i=0; i<sz; i++){
      //std::cout << "el num nodes: " << elements[i]->NumNodes() << " " <<  elements[i]->getOrder() << " " << elements[i]->getElementType() << std::endl;
      temp.emplace_back(elements[i],idealElements.at(elements[i]));
    }

    double merit0 = 0;
    for(int i = 0; i < elements.size(); i++){
      merit0+= all_merits[elements[i]];
    }

    //std::cout << "nodes: " << std::endl;
    //std::cout << nodes.at(nd->first)->xyzvec3() << std::endl;

    assert(elements[0] != elements[1]);

    SubmeshMeritEvaluator evaluator(nodes.at(nd->first).get(),
				    temp,
				    optel_manager, all_merits, factor);


    LBFGSOptimizer optimizer(evaluator);

    optimizer.Optimize();
      


    //std::cout << "After optimize" << std::endl;

    const std::vector<double>& merits = evaluator.getMerits();
    for(int i=0; i<merits.size(); i++){
      all_merits[elements[i]] = merits[i];
    }

    double newmerit = 0;
    for(int i = 0; i < elements.size(); i++){
      newmerit += merits[i];
    }
    //assert(newmerit < merit0);

    }
  }


  double glob_merit = 0.0;
  int cnt=0;
  for(auto el = all_merits.begin(); el != all_merits.end(); ++el){
    if(el->first->getDim() == el_dim_to_opt){
      double fac = 1.0;
      //if(el->first->IsBL()) fac = factor;
      glob_merit+= pow(fac*(el->second-1.0),2);
      //std::cout << el->second << std::endl;
      cnt++;
    }
  }
  return glob_merit/cnt;

  //return activenodes.size();
}
int SubmeshOptimizer::Optimize(double threshold, std::vector<bool> to_opt,
			       std::set<unsigned short int> geo_faces_to_opt_t,
			       std::unordered_set<int> clamped_nodes_t){

  geo_faces_to_opt = geo_faces_to_opt_t;

  //only_opt_node_with_bctag = only_opt_node_with_bctag_t;

  dims_to_opt = to_opt;

  clamped_nodes = clamped_nodes_t;

  using std::cout;
  using std::endl;



  /*
  const node_map& nodes = mesh.getNodes();
  for(auto nd = nodes.begin(); nd != nodes.end(); ++nd){
    cout << nd->second->getType() << endl;
  }
  */
  cout << "in  optimize one by one" << endl;

  //RandomizeNodes(0.01);
  //return 1;

  arma::wall_clock timer;
  timer.tic();
  arma::wall_clock inner_timer;

 
  cout << "initial min Mesh quality: " << minMeshQuality() << endl;

  int max_it = 20;
  //int Nactive = 0;
  //int Nactive_old;
  double Merit = 1.0e100;
  double Merit_old;
  for(int i=0; i<max_it; i++){
    inner_timer.tic();
    //Nactive_old = Nactive;
    Merit_old = Merit;
    //Merit = OptimizeGlobal(0,0,-1,1);
    Merit = OptimizeOneByOne(threshold,1,-1,10);

    //OptimizeGlobal(threshold,2,-1,1.0e0);
    cout << "Merit" << Merit << endl;
    double diff = std::abs(Merit_old-Merit);
    //double diff = 1.0;
    cout << "diff: " << diff << endl;
    if(diff < 1.0e-3){
      //if(diff < 1.0e-2){
      //cout << "Breaking" << endl;
      break;
    }
    //if(Merit != Merit) return 0;
 
    
    //OptimizeOneByOne(threshold,0,-1,std::min(double(i)/2.0+1,4.0));
    //OptimizeGlobal(threshold,0,-1,4.0);
    //OptimizeOneByOne(threshold,0,-1,8.0);
    //OptimizeGlobal(threshold,0,-1,std::min(double(i)/2.0+1,4.0));
    cout << "1x1 time: " << inner_timer.toc() << endl;

  }
  UpdateMeshMerits();

  //OptimizeOneByOne(threshold,0,3,1.0e4);
  inner_timer.tic();
  //SnapToBLExtent();
  if(dims_to_opt[mesh.MeshDimension()]){
    //SnapToBLExtent();
  }
  cout << "min Mesh quality: " << minMeshQuality() << endl;

  cout << "Snap time: " << inner_timer.toc() << endl;
  cout << "BL optimiztion time: " << timer.toc() << endl; 

  
}

/*
int SubmeshOptimizer::SnapToBLExtent(){

  std::unordered_set<MEl*> watch;
  std::map<gind,gind> activenodes;

  const node_map& nodes = mesh.getNodes();
  std::map<MEl*,double>& all_merits = mesh.getAllMeritsNC();
  const ideal_map& idealElements = mesh.getIdealElements();
  const element_set& elements = mesh.getElements();

  std::map<gind, std::vector<MEl*> > node2elements;
  std::unordered_set<gind> bl_edge_nodes;

  const double min_qual0 = minMeshQuality();
  cout << "min qual 0: " << min_qual0 << endl;


  for(auto nd = newnode_map.begin(); nd != newnode_map.end(); ++nd){
    bl_edge_nodes.insert(nd->second);
  }

 
  watch.clear();
  for(auto el = elements.begin(); el != elements.end(); ++el){
    MEl* currel = el->get();
    if(dims_to_opt[currel->getDim()]){
      const int ncn = currel->numCornerNodes();
      const gind* cn = currel->getCornerNodes();
      for(int i=0; i<ncn; i++){
	auto it = bl_edge_nodes.find(cn[i]);
	if(it != bl_edge_nodes.end()){
	  if(1){
	    //if(dims_to_opt[nodes.at(*it)->getType()]){
	    watch.insert(currel);
	    node2elements[*it].push_back(currel);
	  }
	}
      }
    }
  }

    std::map<gind,arma::vec3> oldpos;
  
    for(auto nd = normal_map.begin(); nd != normal_map.end(); ++nd){
	const gind sn = nd->first;
	const gind bln = newnode_map.at(nd->first);
	MNode* BLNode = nodes.at(bln).get();
	MNode* SurfNode = nodes.at(sn).get();
	if(1){
	  //if(dims_to_opt[BLNode->getType()]){
	  arma::vec3 pos = SurfNode->getXYZ();
	  oldpos[bln] = BLNode->getXYZ();
	  pos-= nd->second.dist*nd->second.normal;
	  BLNode->setXYZ(pos.memptr());
	}
    }

    arma::wall_clock timer;
    timer.tic();
    UpdateMeshMerits(watch);
    cout << "Update merit time: " << timer.toc() << endl;


    while(minMeshQuality() < 0.5*min_qual0){
      //for(auto nd = normal_map.begin(); nd != normal_map.end(); ++nd){
      //const gind bln = newnode_map.at(nd->first);
      for(auto nd = node2elements.begin(); nd != node2elements.end(); ++nd){
	const std::vector<MEl*>& els = nd->second;
	
	bool updated = false;
	for(auto el = els.begin(); el != els.end(); ++el){
	  if(all_merits[*el] > 2.0/min_qual0){
	    MNode* BLNode = nodes.at(nd->first).get();
	    BLNode->setXYZ(oldpos[nd->first].memptr());
	    updated = true;
	    break;

	  }
	}
	if(updated){
	  arma::mat gradMerit;
	  for(auto el = els.begin(); el != els.end(); ++el){
	    std::unique_ptr<OptEl> optel = 
	      optel_manager.CreateOptEl(*el,idealElements.at(*el));
	    double detS;
	    optel->computeGradMerit(gradMerit,all_merits[*el],detS);
	  }

	}

      }
      }
    cout << "Min quality: " << minMeshQuality() << endl;
    
      cout << "after moving back" << endl;
 
}
*/
