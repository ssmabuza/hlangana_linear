
#ifndef MassMatrix_HPP
#define MassMatrix_HPP

#include <Teuchos_RCP.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "ParGrid2D.hpp"

class MassMatrix
{

public:

  //! Default Ctor
  MassMatrix(const Teuchos::RCP<ParGrid2D>& par_grid, Teuchos::RCP<Tpetra::CrsMatrix<double>>& mass_mat);

  //! Dtor
  virtual ~MassMatrix() {}
  
  //! Evaluate
  void evaluate();

  //! Extract results view
  const Teuchos::RCP<Tpetra::CrsMatrix<double>>& getvalue() const;

private:

  //! global mesh
  const Teuchos::RCP<ParGrid2D> m_grid;

  //! mass matrix pointer
  Teuchos::RCP<Tpetra::CrsMatrix<double>> m_mat;
  
  //! flag true after computing
  bool m_is_computed;

};

#endif
