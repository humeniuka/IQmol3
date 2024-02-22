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

#include "Data/AtomicProperty.h"
#include "Data/Energy.h"
#include "Data/GeometryList.h"
#include "MoldenParser.h"
#include "TextStream.h"
#include "Util/Constants.h"

#include <QtDebug>
#include <QRegularExpression>

namespace IQmol {
namespace Parser {

bool Molden::parse(TextStream& textStream)
{
   Data::Geometry* currentGeometry(0);
   // Normal modes are given relative to a reference geometry.
   Data::Geometry* referenceGeometry(0);
   Data::VibrationalModeList vibrationalModeList;

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
      line = textStream.previousLine();

      /// Read atoms and coordinates.
      if (line.contains("[Atoms]")) {
         currentGeometry = readAtoms(textStream);
         if (currentGeometry) {
            m_dataBank.append(currentGeometry);
         }
      }
      /// Read frequencies and corresponding normal coordinates.
      else if (line.contains("[FREQ]")) {
         readVibrationalFrequencies(textStream, vibrationalModeList);
      } else if (line.contains("[FR-COORD]")) {
         referenceGeometry = readReferenceCoordinates(textStream);
         if (referenceGeometry) {
            m_dataBank.append(referenceGeometry);
         }
      } else if (line.contains("[FR-NORM-COORD]")) {
         readNormalModes(textStream, vibrationalModeList);
      } else if (line.contains("[INT]")) {
         readVibrationalIntensities(textStream, vibrationalModeList);
      } else if (line.contains("[MASSES]")) {
         readAtomMasses(textStream, referenceGeometry ? referenceGeometry: currentGeometry);
      } else {
         textStream.nextLine();
      }
   }

   // Build frequencies object with information from [FREQ], [INT], [MASSES] and [FR-NORM-COORD] sections.
   if (vibrationalModeList.size() > 0) {
      Data::Frequencies* frequencies = new Data::Frequencies;
      Data::VibrationalModeList::const_iterator iter;
      for (iter = vibrationalModeList.begin(); iter != vibrationalModeList.end(); ++iter) {
         Data::VibrationalMode *vibration = *iter;
         // mass-weighted normal mode displacement
         Data::VibrationalMode *mw_vibration = new Data::VibrationalMode;
         // copy data
         mw_vibration->setFrequency(vibration->frequency());
         mw_vibration->setIntensity(vibration->intensity());
         mw_vibration->setRamanIntensity(vibration->ramanIntensity());

         // BAGEL and QChem do not save the normal mode displacements without the mass-weighting.
         // Since only the mass-weighted normal modes are orthogonal, we apply the mass-weighting again.
         Data::Geometry *geometry = referenceGeometry;
         if (geometry) {
            for (int atom = 0; atom < geometry->nAtoms(); atom++)
            {
               qglviewer::Vec const & vec = vibration->eigenvector()[atom];
               double mass = geometry->getAtomicProperty<Data::Mass>(atom).value();
               double sqrt_mass = sqrt(mass);

               // [dx,dy,dz] --> sqrt(mass) * [dx,dy,dz]
               qglviewer::Vec mw_vec(
                  sqrt_mass * vec.x,
                  sqrt_mass * vec.y,
                  sqrt_mass * vec.z);
               mw_vibration->appendDirectionVector(mw_vec);
            }
            delete vibration;
            vibration = mw_vibration;
         }
         frequencies->append(vibration);
      }
      m_dataBank.append(frequencies);
   }

   return m_errors.isEmpty();
}

bool Molden::lineContainsSectionHeader(QString& line)
{
   /// Check if the line contains a section header,
   /// such as [Atoms], [GTO], [MO] etc.
   const QRegularExpression section_header_rx("\\[[\\-\\w]+\\]");
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

Data::Geometry* Molden::readReferenceCoordinates(TextStream& textStream) {
   Data::Geometry* geometry(new Data::Geometry);
   /// The reference geometry for the normal mode displacements
   /// is given in the format
   ///    [FR-COORD]
   ///    atom_1_element_string x y z
   ///    ...
   ///    atom_n_element_string x y z

   QStringList tokens;
   QString line;

   // Initialize these in case the parser is used more than once.
   m_errors.clear();

   int atomic_number;
   double x, y, z;
   int number_of_atoms = 0;

   while (!textStream.atEnd()) {
      line = textStream.nextNonEmptyLine();
      if (lineContainsSectionHeader(line)) {
         // Begin of next section [...]
         break;
      }

      tokens = TextStream::tokenize(line);
      bool allOk(tokens.size() == 4), ok;

      if (allOk) {
         atomic_number = Data::Atom::atomicNumber(tokens[0]);
         allOk = allOk && (atomic_number != 0);
         x = tokens[1].toDouble(&ok); allOk = allOk && ok;
         y = tokens[2].toDouble(&ok); allOk = allOk && ok;
         z = tokens[3].toDouble(&ok); allOk = allOk && ok;
      }

      if (allOk) {
         geometry->append(atomic_number, qglviewer::Vec(x,y,z));
         number_of_atoms += 1;
      } else {
         QString message("Invalid format line ");
         message += QString::number(textStream.lineNumber());
         message += ":\nExpected: element_string   x  y  z";
         m_errors.append(message);
         break;
      }
   }

   if (number_of_atoms == 0) {
      delete geometry;
      geometry = 0;
   } else {
      // Coordinates in [FR-COORD] are given in Bohr.
      geometry->scaleCoordinates(Constants::BohrToAngstrom);
   }

   return geometry;
}

void Molden::readVibrationalFrequencies(
      TextStream& textStream,
      Data::VibrationalModeList& vibrationalModeList
   ) {
   /// Vibrational frequencies are given on cm^-1. The format is
   ///
   ///   [FREQ]
   ///   frequency_1
   ///   ...
   ///   frequency_n
   ///
   QStringList tokens;
   QString line;

   // If an empy list of vibrations is given, we should build the list,
   // otherwise there is already a list where the data can be set.
   bool createVibrations = (vibrationalModeList.size() == 0);

   // count the number of modes being read.
   int number_of_modes = 0;

   while (!textStream.atEnd()) {
      line = textStream.nextNonEmptyLine();
      if (lineContainsSectionHeader(line)) {
         // Begin of next section [...]
         break;
      }

      double frequency;

      tokens = TextStream::tokenize(line);
      bool allOk(tokens.size() == 1), ok;

      if (allOk) {
         frequency = tokens[0].toDouble(&ok); allOk = allOk && ok;
      }

      if (allOk) {
         if (createVibrations) {
            // Create a new vibration.
            Data::VibrationalMode* vibration = new Data::VibrationalMode(frequency);
            vibrationalModeList.push_back(vibration);
         } else {
            // List of vibrational modes already exists.
            if (number_of_modes >= vibrationalModeList.size()) {
               QString message("Too many vibrational modes in [FREQ] section, line ");
               message += QString::number(textStream.lineNumber());
               m_errors.append(message);
               break;
            }
            Data::VibrationalMode* vibration = vibrationalModeList[number_of_modes];
            vibration->setFrequency(frequency);
         }
         number_of_modes += 1;
      } else {
         QString message("Invalid format line ");
         message += QString::number(textStream.lineNumber());
         message += ":\nExpected one frequency (in cm-1) per line.";
         m_errors.append(message);
         break;
      }
   }
}

void Molden::readNormalModes(
      TextStream& textStream,
      Data::VibrationalModeList& vibrationalModeList
   ) {
   /// Normal mode displacement vectors are in atomic units (Bohr) and
   /// follow the format
   ///      [FR-NORM-COORD]
   ///   vibration vibration_number_1
   ///   atom_1_dx atom_1_dy atom_1_dz
   ///   ...
   ///   atom_n_dx atom_n_dy atom_n_dz
   ///   ....
   ///   vibration vibration_number_N
   ///   atom_1_dx atom_1_dy atom_1_dz
   ///   ...
   ///   atom_n_dx atom_n_dy atom_n_dz

   QStringList tokens;
   QString line;

   // If an empy list of vibrations is given, we should build the list,
   // otherwise there is already a list where the data can be set.
   bool createVibrations = (vibrationalModeList.size() == 0);

   Data::VibrationalMode *vibration(0);
   int atom_index = 0;
   int number_of_modes = 0;

   while (!textStream.atEnd()) {
      line = textStream.nextNonEmptyLine();
      if (lineContainsSectionHeader(line)) {
         // Begin of next section [...]
         break;
      }

      if (line.contains("vibration")) {
         if (createVibrations) {
            // Create a new vibration.
            vibration = new Data::VibrationalMode();
            vibrationalModeList.push_back(vibration);
         } else {
            // List of vibrational modes already exists.
            if (number_of_modes >= vibrationalModeList.size()) {
               QString message("Too many vibrational modes in [FR-NORM-COORD] section, line ");
               message += QString::number(textStream.lineNumber());
               m_errors.append(message);
               break;
            }
            vibration = vibrationalModeList[number_of_modes];
         }
         // Start a new vibrational mode.
         number_of_modes += 1;
         // Start counting atoms from zero again for the new mode.
         atom_index = 0;
         continue;
      }

      // Read displacement vectors for the current vibrational mode.
      if (vibration) {
         tokens = TextStream::tokenize(line);
         bool allOk(tokens.size() == 3), ok;
         double dx, dy, dz;
         // displacement vectors are in Bohr.
         double s = Constants::BohrToAngstrom;

         if (allOk) {
            dx = tokens[0].toDouble(&ok) * s; allOk = allOk && ok;
            dy = tokens[1].toDouble(&ok) * s; allOk = allOk && ok;
            dz = tokens[2].toDouble(&ok) * s; allOk = allOk && ok;
         }

         if (allOk) {
            vibration->appendDirectionVector(qglviewer::Vec(dx, dy, dz));
            atom_index += 1;
         } else {
            QString message("Invalid format line ");
            message += QString::number(textStream.lineNumber());
            message += ":\nExpected 3 floats:  dx  dy  dz";
            m_errors.append(message);
            break;
         }
      }
   }
}

void Molden::readVibrationalIntensities(
      TextStream& textStream,
      Data::VibrationalModeList& vibrationalModeList
   ) {
   /// Vibrational intensities are given in the format
   ///
   ///   [INT]
   ///   ir_intensity_1 [raman_intensity_1]
   ///   ...
   ///   ir_intensity_n [raman_intensity_n]
   ///

   QStringList tokens;
   QString line;

   // If an empy list of vibrations is given, we should build the list,
   // otherwise there is already a list where the data can be set.
   bool createVibrations = (vibrationalModeList.size() == 0);

   // count the number of modes being read.
   int number_of_modes = 0;

   while (!textStream.atEnd()) {
      line = textStream.nextNonEmptyLine();
      if (lineContainsSectionHeader(line)) {
         // Begin of next section [...]
         break;
      }

      double IRintensity = -1.0;
      double ramanIntensity = -1.0;

      tokens = TextStream::tokenize(line);
      bool allOk(tokens.size() <= 2), ok;
      bool hasRaman = (tokens.size() == 2);

      if (allOk) {
         IRintensity = tokens[0].toDouble(&ok); allOk = allOk && ok;
         if (hasRaman) {
            ramanIntensity = tokens[1].toDouble(&ok); allOk = allOk && ok;
         }
      }

      if (allOk) {
         if (createVibrations) {
            // Create a new vibration.
            Data::VibrationalMode* vibration = new Data::VibrationalMode();
            vibration->setIntensity(IRintensity);
            if (hasRaman) {
               vibration->setRamanIntensity(ramanIntensity);
            }
            vibrationalModeList.push_back(vibration);
         } else {
            // List of vibrational modes already exists.
            if (number_of_modes >= vibrationalModeList.size()) {
               QString message("Too many vibrational modes in [INT] section, line ");
               message += QString::number(textStream.lineNumber());
               m_errors.append(message);
               break;
            }
            Data::VibrationalMode* vibration = vibrationalModeList[number_of_modes];
            vibration->setIntensity(IRintensity);
            if (hasRaman) {
               vibration->setRamanIntensity(ramanIntensity);
            }
         }
         number_of_modes += 1;
      } else {
         QString message("Invalid format line ");
         message += QString::number(textStream.lineNumber());
         message += ":\nExpected one or two floats: ir_intensity [raman_intensity].";
         m_errors.append(message);
         break;
      }
   }
}

void Molden::readAtomMasses(
      TextStream& textStream,
      Data::Geometry *geometry
   ) {
   /// Atomic masses are given in the format
   ///
   ///   [MASSES]
   ///   O      15.9949146400
   ///   C      12.0000000000
   ///   N      14.0030740080
   ///   C      12.0000000000
   ///   O      15.9949146400
   ///   N      14.0030740080
   ///   C      12.0000000000
   ///   C      12.0000000000
   ///
   /// The [MASSES] field is not part of the original Molden format specification,
   /// but it is used by TeraChem.
   QStringList tokens;
   QString line;
   QList<double> masses;

   while (!textStream.atEnd()) {
      line = textStream.nextNonEmptyLine();
      if (lineContainsSectionHeader(line)) {
         // Begin of next section [...]
         break;
      }

      tokens = TextStream::tokenize(line);
      bool allOk(tokens.size() == 2), ok;

      if (allOk) {
         double mass = tokens[1].toDouble(&ok); allOk = allOk && ok;
         masses.append(mass);
      } else {
         QString message("Invalid format line ");
         message += QString::number(textStream.lineNumber());
         message += ":\nExpected:  element  mass";
         m_errors.append(message);
         break;
      }
   }
   if (!geometry) {
      QString message("The [MASSES] section has to come after the geometry definition.");
      m_errors.append(message);
      return;
   }

   geometry->setAtomicProperty<Data::Mass>(masses);
}

} } // end namespace IQmol::Parser
