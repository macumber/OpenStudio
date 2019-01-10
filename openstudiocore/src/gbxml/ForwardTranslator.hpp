/***********************************************************************************************************************
*  OpenStudio(R), Copyright (c) 2008-2019, Alliance for Sustainable Energy, LLC, and other contributors. All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
*  following conditions are met:
*
*  (1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following
*  disclaimer.
*
*  (2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
*  disclaimer in the documentation and/or other materials provided with the distribution.
*
*  (3) Neither the name of the copyright holder nor the names of any contributors may be used to endorse or promote products
*  derived from this software without specific prior written permission from the respective party.
*
*  (4) Other than as required in clauses (1) and (2), distributions in any form of modifications or other derivative works
*  may not use the "OpenStudio" trademark, "OS", "os", or any other confusingly similar designation without specific prior
*  written permission from Alliance for Sustainable Energy, LLC.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER(S) AND ANY CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
*  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
*  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER(S), ANY CONTRIBUTORS, THE UNITED STATES GOVERNMENT, OR THE UNITED
*  STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
*  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
*  USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
*  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
*  ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
***********************************************************************************************************************/

#ifndef GBXML_FORWARDTRANSLATOR_HPP
#define GBXML_FORWARDTRANSLATOR_HPP

#include "gbXMLAPI.hpp"

#include "../utilities/core/Path.hpp"
#include "../utilities/core/Optional.hpp"
#include "../utilities/core/Logger.hpp"
#include "../utilities/core/StringStreamLogSink.hpp"

#include "../model/ModelObject.hpp"

#include <map>

class QDomDocument;
class QDomElement;
class QDomNodeList;

namespace pugi {
  class xml_node;
}

namespace openstudio {

  class ProgressBar;
  class Transformation;

namespace model {
  class Model;
  class ModelObject;
  class Material;
  class ConstructionBase;
  class Facility;
  class Building;
  class BuildingStory;
  class ThermalZone;
  class Space;
  class ShadingSurfaceGroup;
  class BuildingStory;
  class Surface;
  class SubSurface;
  class ShadingSurface;
}

namespace gbxml {

  class GBXML_API ForwardTranslator {
  public:

    ForwardTranslator();

    virtual ~ForwardTranslator();

    bool modelToGbXML(const openstudio::model::Model& model, const openstudio::path& path, ProgressBar* progressBar = nullptr);

      /** Get warning messages generated by the last translation. */
    std::vector<LogMessage> warnings() const;

    /** Get error messages generated by the last translation. */
    std::vector<LogMessage> errors() const;

  private:

    QString escapeName(const std::string& name);
    std::string escapeNameS(const std::string& name);

    // listed in translation order
    boost::optional<QDomDocument> translateModel(const openstudio::model::Model& model);
    boost::optional<QDomElement> translateFacility(const openstudio::model::Facility& facility, QDomDocument& doc);
    boost::optional<QDomElement> translateBuilding(const openstudio::model::Building& building, QDomDocument& doc);
    boost::optional<QDomElement> translateSpace(const openstudio::model::Space& space, QDomDocument& doc);
    boost::optional<QDomElement> translateShadingSurfaceGroup(const openstudio::model::ShadingSurfaceGroup& shadingSurfaceGroup, QDomDocument& doc);
    boost::optional<QDomElement> translateBuildingStory(const openstudio::model::BuildingStory& story, QDomDocument& doc);
    boost::optional<QDomElement> translateSurface(const openstudio::model::Surface& surface, QDomDocument& doc);
    boost::optional<QDomElement> translateSubSurface(const openstudio::model::SubSurface& subSurface, const openstudio::Transformation& transformation, QDomDocument& doc);
    boost::optional<QDomElement> translateShadingSurface(const openstudio::model::ShadingSurface& shadingSurface, QDomDocument& doc);
    boost::optional<QDomElement> translateThermalZone(const openstudio::model::ThermalZone& thermalZone, QDomDocument& doc);
    boost::optional<QDomElement> translateLayer(const openstudio::model::Material& material, QDomDocument& doc);
    boost::optional<QDomElement> translateMaterial(const openstudio::model::Material& material, QDomDocument& doc);
    boost::optional<QDomElement> translateConstructionBase(const openstudio::model::ConstructionBase& constructionBase, QDomDocument& doc);
    boost::optional<QDomElement> translateCADObjectId(const openstudio::model::ModelObject& modelObject, QDomElement& parentElement, QDomDocument& doc);

    boost::optional<pugi::xml_node> translateLayer(const openstudio::model::Material& material, pugi::xml_node& root);
    boost::optional<pugi::xml_node> translateMaterial(const openstudio::model::Material& material, pugi::xml_node& root);
    boost::optional<pugi::xml_node> translateConstructionBase(const openstudio::model::ConstructionBase& constructionBase, pugi::xml_node& root);
    boost::optional<pugi::xml_node> translateCADObjectId(const openstudio::model::ModelObject& modelObject, pugi::xml_node& parentElement);

    std::map<openstudio::Handle, QDomElement> m_translatedObjects;
    std::map<openstudio::Handle, pugi::xml_node> m_translatedObjectsS;

    std::set<openstudio::model::Material, openstudio::IdfObjectImplLess> m_materials;

    StringStreamLogSink m_logSink;

    ProgressBar* m_progressBar;

    REGISTER_LOGGER("openstudio.gbxml.ForwardTranslator");
  };

} // gbxml
} // openstudio

#endif // GBXML_FORWARDTRANSLATOR_HPP
