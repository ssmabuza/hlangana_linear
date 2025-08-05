// @HEADER
// ****************************************************************************
//                Hlangana: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hlangana_TrilinosSolverConfig_HPP__
#define __Hlangana_TrilinosSolverConfig_HPP__

#include "Hlangana_OptionHandler.hpp"

namespace hlangana
{

  class TrilinosSolverConfig
      : public Optionable
  {
  public:
    TrilinosSolverConfig(const std::shared_ptr<OptionHandler> &optionHandler)
        : Optionable(optionHandler),
          mOptionHandler(optionHandler)
    {
    }

    ~TrilinosSolverConfig() {}

  private:
    std::shared_ptr<OptionHandler> mOptionHandler;
  };

}
// namespace hlangana

#endif /** __Hlangana_TrilinosSolverConfig_HPP__ */