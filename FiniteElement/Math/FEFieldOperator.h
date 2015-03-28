#pragma once
#include "Field.h"
#include <armadillo>

// LocalState local_state = assembler.GetLocalState(el,sf);

// LocalU u = local_state.u();
// LocalQ q = local_state.q();


// Equation R1 = q + Grad(u);
// Equation R2 = -Div(q) - MyF();
// Assembler.assemble(R1,R2);


// ElementBlock A = Laplacian(u);
// ElementBlock B = Div(Sigma(u,q));
// ElementBlock C = Div(Grad(u));
// ElementBlock E = Grad(Div(u));

namespace FEFieldOperator{
  arma::mat** Laplacian(const arma::mat& u){}



}


