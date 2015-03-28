#include "GlobalMeshOptimizer.h"
#include "LocalMeshMeritEvaluator.h"
#include "NewtonOptimizer.h"

int GlobalMeshOptimizer::Optimize(double threshold, std::vector<bool> to_opt){
  using namespace std;
  /*
  dims_to_opt = to_opt;

  std::unordered_set<MEl*> watch;
  std::map<gind,gind> activenodes;
  const ideal_map& idealElements = mesh.getIdealElements();

  cout << "Threshold: " << threshold << endl;
  FindActiveElements(watch,activenodes,threshold,3,-1);
  cout << "watch size: " << watch.size() << endl;
  cout << "ideal size: " << idealElements.size() << endl;
  cout << "Active nodes size: " << activenodes.size() << endl;

  std::vector<MEl*> temp(watch.size());
  int cnt = 0;
  for(auto el = watch.begin(); el != watch.end(); ++el, ++cnt){
    temp[cnt] = *el;
  }
  

  LocalMeshMeritEvaluator evaluator(activenodes,temp,
				       idealElements,optel_manager,1.0);
  */

  /*

  NewtonOptimizer optimizer(evaluator);

  arma::wall_clock timer;
  timer.tic();
  optimizer.Optimize();
  cout << "Newton optimization time: " << timer.toc() << endl;

  std::map<MEl*,double>& all_merits = mesh.getAllMeritsNC();
  {
    int cnt = 0;
    const std::vector<double>& merits = evaluator.getMerits();
    for(auto el = watch.begin(); el != watch.end(); ++el, ++cnt){
      all_merits[*el] = merits[cnt];
    }
  }
  */

}
