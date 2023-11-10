#ifndef IQMOL_PARSER_MOLDEN_H
#define IQMOL_PARSER_MOLDEN_H
/*******************************************************************************
       
  Copyright (C) 2022 Andrew Gilbert
           
  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.
       
  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  IQmol is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.
      
  You should have received a copy of the GNU General Public License along
  with IQmol.  If not, see <http://www.gnu.org/licenses/>.  
   
********************************************************************************/

#include "Parser.h"
#include "Data/Geometry.h"
#include "Data/Frequencies.h"
#include "Data/VibrationalMode.h"

class QString;

namespace IQmol {
namespace Parser {

   /// Parser for files in Molden format, which is described in
   ///   https://www.theochem.ru.nl/molden/molden_format.html
   /// The Molden format is used by a number of quantum-chemistry packages to store
   /// geometries, optimization trajectories, molecular orbitals and vibrational modes.
   class Molden : public Base {

      public:

         bool parse(TextStream&);

      private:

         bool lineContainsSectionHeader(QString& line);

         Data::Geometry* readAtoms(TextStream&);
         void readAtomMasses(TextStream&, Data::Geometry*);

         Data::Geometry* readReferenceCoordinates(TextStream&);
         void readVibrationalFrequencies(TextStream&, Data::VibrationalModeList&);
         void readNormalModes(TextStream&, Data::VibrationalModeList&);
         void readVibrationalIntensities(TextStream&, Data::VibrationalModeList&);

   };

} } // end namespace IQmol::Parser

#endif
