#include "BoundaryLayer.h"
#include "SubmeshOptimizer.h"
//#include "OptElManager.h"
#include "MeshContainer.h"
#include "GlobalMeshOptimizer.h"


//#include "ShapeFunctionMatrices.h"
//#include "MeshContainer.h"
//#include "NodeIndexer.h"
//#include "ActiveMEl.h"
//#include "OptEl.h"
//#include "ShapeFunctionMatrices.h"
//#include <lbfgs.h>


int BoundaryLayerGenerator::OptimizeBL(){

  using namespace std;

  //ShapeFunctionMatricesFactory sf_factory;

  
  OptElManager optel_manager(mesh.getNodesNC(),geometry,
			     sf_factory,index_factory);



  //SubmeshOptimizer mesh_optimizer(mesh,optel_manager,ExtrudedIdeal,
  //				  normal_map, newnode_map);

  
  //lobalMeshOptimizer global_optimizer(mesh,optel_manager);


  //LocalMeshMeritEvaluator merit_evaluator(
  //LBFGSOptimizer lbfgs_optimizer(merit_evaluator);
  cout << "after mesh_optimizer constructor" << endl;


  /*
  std::vector<bool> to_opt(4,false);
  for(int i = 1; i <= mesh.MeshDimension(); i++){
    to_opt[i] = true;
  }
  mesh_optimizer.Optimize(1.2,to_opt);
  */
  /*
  for(int i = 1; i <= mesh.MeshDimension(); i++){
    std::vector<bool> to_opt(4,false);
    to_opt[i] = true;
    mesh_optimizer.Optimize(1.2,to_opt);
  }
  */


  //SubmeshOptimizer mesh_optimizer(mesh,optel_manager,ideal_elements,
  //				  2); 
//mesh_optimizer.Optimize(1.2,{0,0,1,0},symmetry_faces);

  SubmeshOptimizer mesh_optimizer(mesh,optel_manager,ideal_elements,
  				  mesh.Dimension()); 

  if(mesh.Dimension() == 3){
    mesh_optimizer.Optimize(1.2,{0,0,1,1},symmetry_faces,clamped_nodes);
  }
  else{
    mesh_optimizer.Optimize(1.2,{0,0,1,0});
  }
  /*
  for(int dim = 3; dim <= mesh.Dimension(); dim++){

    std::vector<bool> to_opt(4,false);
    to_opt[dim] = true;
    to_opt[2] = true;
    //to_opt[2] = true;
    //to_opt[2] = true;

    
    mesh_optimizer.Optimize(1.2,to_opt,mesh.SymmetrySurface());
 
    //to_opt[3] = false;
    //to_opt[2] = true;
    //mesh_optimizer.Optimize(1.2,to_opt); 

    //to_opt[2] = false;
    //to_opt[3] = true;
    //mesh_optimizer.Optimize(1.2,to_opt);
  }
  */
  //mesh_optimizer.Optimize(0.9,{false,false,true,false});
  //mesh_optimizer.Optimize(1.2,{false,false,false,true});
  //std::vector<bool> to_opt(4,false);
  //to_opt[mesh.Dimension()] = true;
  //mesh_optimizer.Optimize(1.2,to_opt);
  //mesh_optimizer.Optimize(1.2,{false,false,true,true});

  //global_optimizer.Optimize(1.2,{false,false,true,true});

  /*
  IncrementalMeshOptimizer optimizer(mesh,optel_manager,ExtrudedIdeal,
  				     normal_map, newnode_map);
  optimizer.Optimize(1.3);
  */

  /*
  double distortion_threshold = 1.0;
  double merit_threshold = pow(distortion_threshold-1.0,2);
  
  MeritEvaluator merit_evaluator(mesh);

  merit_evaluator.Initialize(4,merit_threshold,-1);
  
 
  cout << "Merit 0: " << merit_evaluator.EvaluateMerit() << endl;

  //return 1;


  ObjectiveFunction objective_function(&merit_evaluator);

  objective_function.Run();
  cout << "optimize time: " << timer.toc() << endl;

  */
}

/*
class ElOptData{
private:
 
public:
 
  ElOptData(const arma::mat& I): ideal(I), merit(1.0/1.0e-14) {}
  const arma::mat& ideal;
  
  double merit;
  //const double getMerit() const{ return merit; }
  //void setMerit(double m) { merit = m; }

};
*/

/*
class GlobalSystem{
private:
  arma::mat Jacobian;
  //arma::umat ij;
  arma::vec val, RHS;

public:
  void Initialize(const int nn);

  void AddToSystem(const int ncn, const gind* cn, arma::vec& gradMerit,
		   arma::mat& HessianMerit, const std::map<gind,gind>& 
		   node_map);

  void Assemble();

  arma::mat Solve();


};

void GlobalSystem::Initialize(const int nn){

  RHS.resize(nn);
  RHS.zeros();
  Jacobian.resize(nn,nn);
  Jacobian.zeros();

  //  ij.resize(0,2);
  //  val.reset();

}

void GlobalSystem::AddToSystem(const int ncn, const gind* cn, 
			       arma::vec& gradMerit,
			       arma::mat& HessianMerit, 
			       const std::map<gind,gind>& node_map){



  for(auto nd = node_map.begin(); nd != node_map.end(); ++nd){
    //cout << nd->first << " " << nd->second << endl;

  }

  //cout << "norm H-H': " << norm(HessianMerit-HessianMerit.t()) << endl;
  // cout << "before cn loop: " << endl;
  for(int i=0; i<ncn; i++){
    const int nd = cn[i];
    if(node_map.find(nd) != node_map.end()){
      const int loc_node = node_map.at(nd);
      // cout << "nd: " << nd << " " << loc_node << endl;
      for(int j=0; j<3; j++){
	RHS(loc_node*3+j)+= gradMerit(i*3+j);
	//RHS(loc_node*3+j)+= gradMerit(j*3+i);
	for(int k=0; k<ncn; k++){
	  const int ndk = cn[k];
	  if(node_map.find(ndk) != node_map.end()){
	    const int loc_node_k = node_map.at(cn[k]);
	    for(int l=0; l<3; l++){
	      Jacobian(loc_node*3+j, loc_node_k*3+l)+= HessianMerit(i*3+j,k*3+l);
	      //Jacobian(loc_node*3+j, loc_node_k*3+l)+= HessianMerit(i+3*j,k+3*l);
	    }
	  }
	}
      }
    }
  }
 
}

void GlobalSystem::Assemble(){

}

arma::mat GlobalSystem::Solve(){
  const double cnd = arma::cond(Jacobian);
  cout << "Condition NO: " << cnd  << endl;
  Jacobian+= arma::eye(Jacobian.n_rows,Jacobian.n_cols);
  //if(cnd > 1.0e10){
    //Jacobian+= std::min(cnd/1.0e12,1.0)*arma::eye(Jacobian.n_rows,Jacobian.n_cols);
  //}
  //cout << "New Cond NO: " << cond(Jacobian) << endl;

  Jacobian.save("Jacobian",arma::raw_ascii);
  //RHS.save("RHS",arma::raw_ascii);
  arma::vec x = arma::solve(Jacobian,RHS);
  return arma::reshape(x,3,x.n_elem/3);
}
*/

/*
const double computeElementMerit(MEl* element,
			node_map& nodes,
			NodeIndexFactory& index_factory,
			ShapeFunctionMatricesFactory& sf_factory,
				 const arma::mat& ideal){



  const ShapeFunctionMatrices* sf = 
    sf_factory.getShapeFunction(element->getElementType(),
				element->getOrder(),0);

  ActiveMEl activeEl(element,index_factory,nodes,sf);


  if(element->getDim() == 2){
    OptEl2D optel(activeEl,ideal);
    return optel.computeMerit();
  }
  else if(element->getDim() == 3){
    OptEl3D optel(activeEl,ideal);
    return optel.computeMerit();
  }

}
*/

/*
double computeMeshMerit(std::vector<std::pair<MEl*,double> >& current_elements,
			node_map& nodes,
			NodeIndexFactory& index_factory,
			ShapeFunctionMatricesFactory& sf_factory,
			const ideal_map& idealElements,
			std::map<MEl*,double>& all_merits){




  //double merit = 0.0;
  for(auto el=current_elements.begin(); el != current_elements.end(); ++el){
    const ShapeFunctionMatrices* sf = 
      sf_factory.getShapeFunction(el->first->getElementType(),
				  el->first->getOrder(),0);

    ActiveMEl activeEl(el->first,index_factory,nodes,sf);


    if(el->first->getDim() == 2){
      //OptEl2D optel(activeEl,el->second);
      //el->second =  optel.computeMerit();
      //merit+= el->second;
    }
    else if(el->first->getDim() == 3){
      OptEl3D optel(activeEl,idealElements.at(el->first));
      el->second = optel.computeMerit();
      //merit+= el->second;
      all_merits[el->first] = el->second;
      //merit+= el->second;
    }
  
  }

  double max_merit = 0.0;
  double merit = 0.0;
  for(auto el = all_merits.begin(); el != all_merits.end(); ++el){
    merit+= el->second;
    max_merit = std::max(el->second,max_merit);
  }


  return merit;

}

void resetActiveNodes(std::map<MEl*,ElOptData>& current_elements,
		      std::map<gind,bool>& active_nodes){

  active_nodes.clear();

  for(auto el = current_elements.begin(); el != current_elements.end(); ++el){
    if(!el->first) cout << "Element is null!" << endl;
    const int ncn = el->first->numCornerNodes();
    const gind* cn = el->first->getCornerNodes();
    for(int i=0; i<ncn; i++){
      active_nodes[cn[i]] = true;
    }
  }
 
}
void modifyCurrentElements(std::map<MEl*,ElOptData>& current_elements,
			   std::vector<std::pair<MEl*,double> >& check_elements,
			   ideal_map& idealElements,
			   const double merit_threshold)
{


  
  for(auto el = check_elements.begin(); el != check_elements.end(); ++el){
    if(el->second < merit_threshold){
      current_elements.erase(el->first);
    }
    else{
      if(current_elements.find(el->first) == current_elements.end()){
	current_elements.emplace(el->first,
				 ElOptData(idealElements.at(el->first)));
      }
    }
  }
  

}

			   
void findActiveElements(std::map<MEl*,ElOptData>& current_elements,
			std::vector<std::pair<MEl*,double> >& check_elements,
			std::vector<bool>& active_elements,
			node_map& nodes,
			NodeIndexFactory& index_factory,
			ideal_map& idealElements,
			const std::vector<MEl*>& element_vector, 
			std::map<gind,bool>& active_nodes,
			const double merit_threshold){

  // Zero-out check elements
  check_elements.clear();


  int cnt = 0;
  for(auto el = element_vector.begin(); el != element_vector.end(); 
      ++el, ++cnt){

    bool has_active = false;
    const int ncn = (*el)->numCornerNodes();
    const gind* cn = (*el)->getCornerNodes();
    for(int i=0; i<ncn; i++){
      if(active_nodes[cn[i]] == true) { has_active = true; break; }
    }
    if(has_active){
      check_elements.push_back(std::make_pair(*el,1.0/1.0e-14));
	    
   
      auto it = idealElements.find(*el);
      if(it == idealElements.end()){
	// Add to ideal map
	ActiveMEl activeEl(*el,index_factory,nodes,NULL);
	idealElements[*el] = activeEl.getNodesMatrix();

      }
	    
    }
  }
  

  for(auto el = current_elements.begin(); el != current_elements.end(); ++el){
    //check_elements.push_back(std::make_pair(el->first,1.0/1.0e-14));
  }



}

*/



/*
int BoundaryLayerGenerator::OptimizeBL(){
 
  arma::wall_clock timer;

  ShapeFunctionMatricesFactory sf_factory;
  NodeIndexFactory index_factory;
  ideal_map& idealElements = mesh.getIdealElementsNC();

  const double distortion_threshold = 1.01;
  const double merit_threshold = pow(distortion_threshold-1.0,2);

  const element_set& elements = mesh.getElements();
  node_map& nodes = mesh.getNodesNC();

  std::vector<MEl*> element_vector(elements.size());
  std::vector<MNode*> node_vector(nodes.size());
  std::vector<bool> active_elements(elements.size());
  std::map<gind,bool> active_nodes;
  std::vector<int> ncn_vector(elements.size());
  std::vector<const gind*> cn_vector(elements.size());
  std::vector<std::pair<MEl*,double> > check_elements;

  std::map<MEl*,ElOptData> current_elements;
  std::map<MEl*,double> all_merits;

  std::map<gind,arma::vec3> grads, pos_old;
  std::map<gind,gind> current_node_map;

  cout << "in optimize BL!" << endl;

  {
    cout << "Ideal size: " << idealElements.size() << endl;
 

    int cnt = 0;
    for(auto el = elements.begin(); el != elements.end(); ++el, ++cnt){
      element_vector[cnt] = el->get(); 
 
      //arma::mat ideal;
      auto it = idealElements.find(el->get());
      if(it == idealElements.end()){
	
	ActiveMEl activeEl(el->get(),index_factory,nodes,NULL);
	arma::mat ideal = activeEl.getNodesMatrix();
	idealElements[el->get()] = ideal;
	current_elements.emplace(el->get(),idealElements[el->get()]);
	check_elements.push_back(std::make_pair(el->get(),1.0/0.0e-14));
	
      }
      else{
	current_elements.emplace(el->get(),it->second);
	check_elements.push_back(std::make_pair(el->get(),1.0/0.0e-14));
      }

    }
  }
  cout << "before reset active nodes" << endl;
  resetActiveNodes(current_elements,active_nodes);
  cout << "after reset active nodes" << endl;
 
  cout << "current element size: " << current_elements.size() << endl;

  for(auto el = current_elements.begin(); el != current_elements.end(); ++el){
    if(el->first->getElementType() == 6 && el->second.ideal.n_cols != 6){
      cout << "ERROR!" << endl;
    }
  }


 
  // Compute initial merit
  double merit0 = computeMeshMerit(check_elements,nodes,index_factory,
				   sf_factory,idealElements, all_merits);

 cout << "current element size 3: " << current_elements.size() << endl;
 cout << "check element size: " << check_elements.size() << endl;
 cout << setprecision(15) << "Merit 0: " << merit0 << endl;


 //GlobalSystem global_system;

  double step = 1.0;
  int goodsteps = 0;
  bool recompute = true;
  double merit=merit0;
  double meritold;
  int allits = 0;
  arma::mat dx;
  cout << "before starting loop" << endl;


  timer.tic();
  while(goodsteps < 1 && step > 1.0e-50 && current_elements.size() > 0){
    
    if(recompute){

      if(goodsteps > 1){
	modifyCurrentElements(current_elements,check_elements,idealElements,
			      merit_threshold);
      }   

      resetActiveNodes(current_elements,active_nodes);

      // Find/erase  active elements
      findActiveElements(current_elements,check_elements,active_elements,
			 nodes,index_factory,idealElements,element_vector,
			 active_nodes,merit_threshold);
   

      // Initialize global system
 


      // Zero-out gradients
      grads.clear();
      current_node_map.clear();
      {
	int cnt2=0;
	for(auto nd = active_nodes.begin(); nd != active_nodes.end(); ++nd){
	  if(nd->second){

	    if(nodes[nd->first]->getType() == 3){
	      grads[nd->first] = {0.0,0.0,0.0};
	      current_node_map[nd->first] = cnt2;
	      cnt2++;
	    }
	  }
	}
	cout << "cnt2 size: " << cnt2 << endl;
	//global_system.Initialize(cnt2*3);
      }
   

      //cout << "before computing gradients" << endl;
      //cout << "current size: " << current_elements.size() << endl;
      // Actually compute the gradients
      timer.tic();
      for(auto el = current_elements.begin(); el != current_elements.end(); ++el){
	const ShapeFunctionMatrices* sf = 
	  sf_factory.getShapeFunction(el->first->getElementType(),
				      el->first->getOrder(),0);
	
	//cout << "in element loop" << endl;
	ActiveMEl activeEl(el->first,index_factory,nodes,sf);

	arma::mat gradMerit;
	arma::mat HessianMerit;
	if(el->first->getDim() == 2){
	  //OptEl2D optel(activeEl,ideal);
	  //gradMerit = optel.computeGradMerit();
	}
	else if(el->first->getDim() == 3){
	  //cout << "b" << endl;
	  //cout << el->first->getElementType() << " " << el->second.ideal.n_cols << endl;
	  OptEl3D optel(activeEl,el->second.ideal);
	  gradMerit = optel.computeGradMerit();
	  HessianMerit = optel.computeHessianMeritFD();
	  //cout << "a" << endl;
	}



	const int ncn = el->first->numCornerNodes();
	const gind* cn = el->first->getCornerNodes();

	for(int i=0; i<ncn; i++){
	  if(current_node_map.find(cn[i]) != current_node_map.end()){
	    if(nodes[cn[i]]->getType() == 3){
	      grads[cn[i]]+= gradMerit.col(i);
	    }
	  }
	}
	arma::vec gradMeritVec = vectorise(gradMerit);
	//global_system.AddToSystem(ncn,cn,gradMeritVec,HessianMerit,
	//							  current_node_map);

      }// End current element loop
      cout << "Hessian time: " << timer.toc() << endl;

      cout << "after element loop" << endl;
      timer.tic();
      //dx = global_system.Solve();
      cout << "solve time: " << timer.toc() << endl;
    }

 

    pos_old.clear();
    for(auto nd = grads.begin(); nd != grads.end(); ++nd){
 
      arma::vec xyz = nodes[nd->first]->getXYZ();
      //cout << "grad: " << endl;
      //cout << nd->second << endl;
      //cout << "Hess grad: " << endl;
      //cout << dx.col(current_node_map.at(nd->first)) << endl;
      pos_old[nd->first] = xyz;
      if(step < 0.0){
	xyz-= step*nd->second;
      }
      else{
	//xyz-= step*dx.col(current_node_map.at(nd->first));
      }
      nodes[nd->first]->setXYZ(xyz.memptr());
    }
    
    // Compute new merit
    meritold = merit;
    merit = computeMeshMerit(check_elements,nodes,index_factory,
				    sf_factory,idealElements, all_merits);

    //arma::vec::fixed<1> temp;
    //temp(0) = merit;
    if(merit < meritold && merit != arma::datum::inf){
      step = std::min(2.0*step,1.0);
      recompute = true;
      cout  << "GOOD STEP " << goodsteps+1 << ". merit: " << 
	merit << " step: " << step << " " << check_elements.size() << " " <<
	current_elements.size() << endl;
      goodsteps++;
    }
    else{
      for(auto nd = pos_old.begin(); nd != pos_old.end(); ++nd){
	nodes[nd->first]->setXYZ(nd->second.memptr());
      }
   
      //cout << "BAD STEP! merit: " << merit << " step: " << step << " " <<
      //	check_elements.size() << endl;

      recompute = false;
      step = 1.0/2.0*step;
      merit = meritold;
 
    }
    allits++;
  }
  cout << "while loop time: " << timer.toc() << endl;



  //std::map<gind,arma::vec3> currentgrad;
  //std::map<gind,arma::vec3> oldpositions;


 
 
 

  cout << "before compute mesh merit: " << endl;
  timer.tic();
  const double merit0 = computeMeshMerit(elements,nodes,index_factory,
					 sf_factory,idealElements);
  cout << "merit0: " << merit0 << endl;
  cout << "merit computation time: " << timer.toc() << endl;
  
  double step = 1.0;
  int goodsteps = 0;
  bool recompute = true;
  double merit=merit0;

  
  
  while(goodsteps < 200 && step > 1.0e-50){
    if(recompute){
      // Find active elements
      for(auto el = elements.begin(); el != elements.end(); ++el){
	auto it = idealElements.find((*el).get());
	if(it == idealElements.end()){ // If not already active
	  const int ncn = (*el)->numCornerNodes();
	  const gind* cn = (*el)->getCornerNodes();
	  bool newactive = false;
	  for(int i=0; i<ncn; i++){
	    //	    if(nodeactive[cn[i]]->second = true){
	    //  newacrive = true;
	    //  break;
	    //}
	    if(currentgrad.find(cn[i]) != currentgrad.end()){
	      newactive = true;
	      break;
	    }
	  }
      
	  if(newactive){
	    const ShapeFunctionMatrices* sf = 
	      sf_factory.getShapeFunction((*el)->getElementType(),
					  (*el)->getOrder(),0);

	    ActiveMEl activeEl((*el).get(),index_factory,nodes,sf);
	    arma::mat ideal = activeEl.getNodesMatrix();
	    OptEl3D optel(activeEl,ideal);
	    if(optel.computeDistortion() > threshold){
	      idealElements[(*el).get()] = ideal;
	      for(int i=0; i<ncn; i++){
		currentgrad[cn[i]] = {0.0,0.0,0.0};
	      }
	    }
	  }
	}
      }
     
      // Zero-out currentgrad
      for(auto nd = currentgrad.begin(); nd != currentgrad.end(); ++nd){
	nd->second.zeros();
      }
  
      for(auto el=elements.begin(); el != elements.end(); ++el){
	const ShapeFunctionMatrices* sf = 
	  sf_factory.getShapeFunction((*el)->getElementType(),
				      (*el)->getOrder(),0);


      
	ActiveMEl activeEl((*el).get(),index_factory,nodes,sf);

	arma::mat ideal;
	auto it = idealElements.find((*el).get());
	if(it != idealElements.end()) ideal = it->second;

	if(!ideal.is_empty()){
	  arma::mat gradMerit;
	  if((*el)->getDim() == 2){
	    OptEl2D optel(activeEl,ideal);
	    gradMerit = optel.computeGradMerit();
	  }
	  else if((*el)->getDim() == 3){
	    OptEl3D optel(activeEl,ideal);
	    gradMerit = optel.computeGradMerit();
	  }
	  if(!gradMerit.is_finite()){
	    //grad_finite = false;
	    //break;
	  }
	  const int ncn = (*el)->numCornerNodes();
	  const gind* cn = (*el)->getCornerNodes();
	  for(int i=0; i<ncn; i++){
	    if(nodes[cn[i]]->getType() == 3){
	      currentgrad[cn[i]]+= gradMerit.col(i);
	    }
	  }

	} // end if
  
      }// End element loop
    }
    for(auto nd = currentgrad.begin(); nd != currentgrad.end(); ++nd){
      //arma::vec3 xyz = nodes[nd->first]->xyzvec3();
      arma::vec3 xyz = nodes[nd->first]->getXYZ();
 
      oldpositions[nd->first] = xyz;
      xyz-= step*nd->second;
      if(!xyz.is_finite()){
	//cout << "xyz not finite!" << endl;
	//cout << nd->second.t() << endl;
      }
      nodes[nd->first]->setXYZ(xyz.memptr());
 
      

    }
    double meritold = merit;
    merit = computeMeshMerit(elements,nodes,index_factory,
					 sf_factory,idealElements);
 
    arma::vec::fixed<1> temp;
    temp(0) = merit;
    if(merit <= meritold && merit != arma::datum::inf){
      //if(merit < meritold && temp.is_finite()){
      //step = std::min(std::max((meritold-merit),2.0)*step,1.0);
      step = std::min(2.0*step,1.0);
      recompute = true;
      cout  << "GOOD STEP! merit: " << merit << " step: " << step << endl;
    
      goodsteps++;
    }
    else{
      for(auto nd = currentgrad.begin(); nd != currentgrad.end(); ++nd){
	arma::vec3 xyz = nodes[nd->first]->getXYZ();
	nodes[nd->first]->setXYZ(oldpositions[nd->first].memptr());



      }
   
      cout << "BAD STEP! merit: " << merit << " step: " << step << endl;
     
      //cout << "merit old: " << meritold << endl;
 
      recompute = false;
      step = 0.5*step;
      merit = meritold;
 
    }
  }
  
  cout << "At end of BL opt!" << endl;
  
}
*/
  /*
  timer.tic();
  for(auto el = idealElements.begin(); el == idealElements.begin(); ++el){
    const ShapeFunctionMatrices* sf = 
      sf_factory.getShapeFunction(el->first->getElementType(),
				  el->first->getOrder(),0);
    
    ActiveMEl activeEl(el->first,index_factory,nodes,sf);
    OptEl3D optel(activeEl,el->second);
    arma::mat HessDeb = optel.computeHessianMeritFD();
    arma::mat HessAnalytical = optel.computeHessianMerit();
    arma::mat gradMerit = optel.computeGradMerit();
    arma::mat gradDebug = optel.debugGradMerit();

    const double distortion = optel.computeDistortion();
    cout << "distortion: " << distortion << endl;
    cout << "debug Merit: " << endl;
    cout << gradDebug << endl;
    cout << "Hess FD: " << endl;
    cout << HessDeb << endl;
    cout << "Hess Analytical: " << endl;
    cout << HessAnalytical << endl;
    cout << "norm difference: " << 
      norm(HessAnalytical-HessDeb,2.0)/norm(HessDeb,2.0) << endl;
    
    //HessAnalytical.save("HessAnalytical",arma::raw_ascii);
  }
  cout << "All element Hessian time: " << timer.toc() << endl;
  */

 /*

      */
/*
double computeMeshMerit(const element_set& elements,
			node_map& nodes,
			NodeIndexFactory& index_factory,
			ShapeFunctionMatricesFactory& sf_factory,
			const ideal_map& idealElements){
 double merit = 0.0;
  double maxmerit = 0.0;

  for(auto el=idealElements.begin(); el != idealElements.end(); ++el){
    const ShapeFunctionMatrices* sf = 
      sf_factory.getShapeFunction(el->first->getElementType(),
				  el->first->getOrder(),0);

    ActiveMEl activeEl(el->first,index_factory,nodes,sf);


    if(el->first->getDim() == 2){
      OptEl2D optel(activeEl,el->second);
      merit+= optel.computeMerit();
    }
    else if(el->first->getDim() == 3){
      OptEl3D optel(activeEl,el->second);
      merit+= optel.computeMerit();
    }
  }

  return merit;

}
*/
