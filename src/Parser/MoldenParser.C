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

#include "Data/Energy.h"
#include "Data/GeometryList.h"
#include "MoldenParser.h"
#include "TextStream.h"
#include "Util/Constants.h"

#include <QRegularExpression>

namespace IQmol {
namespace Parser {

bool Molden::parse(TextStream& textStream)
{
   Data::Geometry* currentGeometry(0);

   QStringList tokens;
   QString line;

   // Check that the file has the correct header.
   if (!textStream.atEnd()) {
      line = textStream.nextLine();
      if (!line.contains("[MoldenFormat]")) {
         // This is not a molden file.
      }
   }

   while (!textStream.atEnd()) {
      line = textStream.nextLine();

      if (line.contains("[Atoms]")) {
         currentGeometry = readAtoms(textStream);
         if (currentGeometry) {
            m_dataBank.append(currentGeometry);
         }
      }
   }

   return m_errors.isEmpty();
}

bool Molden::lineContainsSectionHeader(QString& line)
{
   /// Check if the line contains a section header,
   /// such as [Atoms], [GTO], [MO] etc.
   const QRegularExpression section_header_rx("\\[\\w+\\]");
   QRegularExpressionMatch match(section_header_rx.match(line));
   return match.hasMatch();
}

Data::Geometry* Molden::readAtoms(TextStream& textStream)
{
   // Read coordinates for electron density and molecular orbitals.
   // The format is:
   //
   //    [Atoms] (Angs|AU)
   //    element_name1 number1 atomic_number1 x1 y1 z1
   //    element_name2 number2 atomic_number2 x2 y2 z2
   //    ...

   // Line containing the [...] section marker.
   QString first_line = textStream.previousLine();

   QStringList tokens;
   QString line;

   Data::Geometry* geometry(new Data::Geometry);
   // Initialize these in case the parser is used more than once.
   m_errors.clear();

   int atomic_number;
   double x, y, z;
   int number_of_atoms(0);

   while (!textStream.atEnd()) {
      line = textStream.nextNonEmptyLine();
      if (lineContainsSectionHeader(line)) {
         // Begin of next section [...]
         break;
      }

      tokens = TextStream::tokenize(line);
      bool allOk(tokens.size() == 6), ok;

      if (allOk) {
         atomic_number = tokens[2].toUInt(&ok);
         x = tokens[3].toDouble(&ok); allOk = allOk && ok;
         y = tokens[4].toDouble(&ok); allOk = allOk && ok;
         z = tokens[5].toDouble(&ok); allOk = allOk && ok;
      }

      if (allOk) {
         geometry->append(atomic_number, qglviewer::Vec(x,y,z));
         number_of_atoms += 1;
      } else {
         QString message("Invalid format line ");
         message += QString::number(textStream.lineNumber());
         message += ":\nExpected: element_name number atomic_number   x  y  z";
         m_errors.append(message);
         break;
      }
   }

   if (number_of_atoms == 0) {
      delete geometry;
      geometry = 0;
   } else {
      // By default it is assumed that units are in Angstrom.
      bool convertFromBohr();
      if (first_line.contains("AU")) {
         geometry->scaleCoordinates(Constants::BohrToAngstrom);
      }
   }

   return geometry;
}

} } // end namespace IQmol::Parser
