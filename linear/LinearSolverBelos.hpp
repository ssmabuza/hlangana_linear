
#ifndef LinearSolverBelos_HPP
#define LinearSolverBelos_HPP


#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

/**
 * \brief A routine for solving a linear system using Belos with Ifpack2 Preconditioners
 * 
 * 
 */
void SolveLinearSystemBelos(const Teuchos::RCP<Tpetra::CrsMatrix<double>>& A, 
                            const Teuchos::RCP<Tpetra::Vector<double>>& rhs, 
                            Teuchos::RCP<Tpetra::Vector<double>>& sol,
                            Teuchos::RCP<Teuchos::ParameterList>& prec_pl,
                            Teuchos::RCP<Teuchos::ParameterList>& solver_pl);






#endif