/***********************************************************************************/
/*  Copyright 2020 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef LIBACDD_H
#define LIBACDD_H

#include <meteoio/IOUtils.h>
#include <meteoio/Config.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/dataClasses/Grid2DObject.h>

#include <string>
#include <vector>

/**
 * @class ACDD
 * @brief This class contains and handles NetCDF Attribute Conventions Dataset Discovery attributes (see 
 * <A href="http://wiki.esipfed.org/index.php?title=Category:Attribute_Conventions_Dataset_Discovery">ACDD</A>).
 * @details A few attributes can get their default value automatically from the data. For the others, some "best efforts" are made in order to keep
 * the whole process as simple as possible. It is however possible to provide some of these attributes from the configuration file, using the
 * following keys:
 *  - NC_CREATOR: the name of the creator of the data set (default: login name);
 *  - NC_EMAIL: the email of the creator;
 *  - NC_KEYWORDS: a list of AGU Index Terms (default: hard-coded list);
 *  - NC_TITLE: a short title for the data set;
 *  - NC_INSTITUTION: the institution providing the data set (default: domain name);
 *  - NC_PROJECT: the scientific project that created the data;
 *  - NC_ID: an identifier for the data set, provided by and unique within its naming authority. Example: DOI, URL, text string, but without white spaces
 *  - NC_NAMING_AUTHORITY: The organization that provides the initial id (see above) for the dataset
 *  - NC_PROCESSING_LEVEL: a textual description of the processing level
 *  - NC_SUMMARY: a paragraph describing the dataset;
 *  - NC_SUMMARY_FILE: a file containing a description of the dataset, it overwrites the value of NC_SUMMARY if present;
 *  - NC_COMMENT: miscellaneous informartion about the dataset;
 *  - NC_ACKNOWLEDGEMENT: acknowledgement for the various types of support for the project that produced this data;
 *  - NC_METADATA_LINK: A URL/DOI that gives more complete metadata;
 *  - NC_LICENSE: describes the license applicable to the dataset;
 *  - NC_PRODUCT_VERSION: Version identifier of the data file or product as assigned by the data creator (default: 1.0).
*/
class ACDD {
	public:
		enum Mode {MERGE, REPLACE, APPEND};
		
		ACDD() : name(), cfg_key(), value() {defaultInit();}
		
		void setUserConfig(const mio::Config& cfg, const std::string& section);
		
		void addAttribute(const std::string& att_name, const std::string& att_value, const std::string& att_cfg_key="", Mode mode=MERGE);
		void addAttribute(const std::string& att_name, const double& att_value, const std::string& att_cfg_key="", const Mode& mode=MERGE);

		void getAttribute(const size_t ii, std::string &att_name, std::string & att_value) const;
		size_t getNrAttributes() const {return name.size();}
		
		void setGeometry(const mio::Grid2DObject& grid, const bool& isLatLon);
		void setGeometry(const std::vector< std::vector<mio::MeteoData> >& vecMeteo, const bool& isLatLon);
		void setGeometry(const mio::Coords& location, const bool& isLatLon);
		
		void setTimeCoverage(const std::vector< std::vector<mio::MeteoData> >& vecMeteo);
		void setTimeCoverage(const std::vector<mio::MeteoData>& vecMeteo);
		
	private:
		void defaultInit();
		size_t find(const std::string& search_name) const;
		
		std::vector<std::string> name, cfg_key, value;
};

#endif
