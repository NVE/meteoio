#ifndef FILTERVALUE_H_INCLUDED
#define FILTERVALUE_H_INCLUDED

#include "FilterBase.h"
#include "FilterBase1Stn.h"
#include "FilterBaseMStn.h"

using namespace std;

/**
 * @class FilterValue
 * @brief Helpers for filtering of a measure's values. 
 *        Either use its static methods, or use it as an enhancer for concrete filters; 
 *        in the latter case, consider using the ready-to-use FilterValue1Stn or FilterValueMStn. 
 *        Design note: this class has been thought as a subclass of FilterBase to be used in a mixin, 
 *                     but with an implementation that avoid problematic multiple inheritance. 
 * @author Florian Hof
 * @date   2009-03-02
 */
class FilterValue {

	public:

		// base definition

		/** 
		* Constructor of the FilterValue helper, based on an abstract filter. 
		*/
		FilterValue(FilterBase* filter);
		virtual ~FilterValue(){};

		// check handling

		virtual void prepareCheck();

		// accessors

		/** 
		* Gets the name (usually a technical abbreviation) of the mesure, 
		* as defined in the config file and recognize by this class. 
		*/
		const string& getMeasureName() const;

		/**
		* Get the measure's value of a meteo data. 
		* The measure to return is determined by the corresponding filter's parameter.
		* @param data   The reference of the MeteoData from which the value will be taken.
		* @return The reference to the requested value; the value is thus editable. 
		*/
		double& getMeasureValue(MeteoData& data) const;
		const double& getMeasureValue(const MeteoData& data) const;

		// helpers

		/**
		* Get the measure's value of a meteo data. 
		* @param data   The reference of the MeteoData from which the value will be taken.
		* @param measurePtr   The reference of which measure of the MeteoData to return (as returned by getMeasurePtr).
		* @return The reference to the requested value; the value is thus editable. 
		*/
		static double& getMeasureValue(MeteoData& data, double MeteoData::* measurePtr);
		static const double& getMeasureValue(const MeteoData& data, const double MeteoData::* measurePtr);

		/** 
		* Gets the pointer of a mesure, given its name. 
		* @param measureName   The name (usually a technical abbreviation) of the mesure, 
		*                      as defined in the config file and recognize by this class. 
		* @return The reference of the measure within a MeteoData, useable by getMeasureValue. 
		*/
		static double MeteoData::* getMeasurePtr(const string& measureName);

	private:

		FilterBase* m_filter;

		// interpreted parameters

		/** The name (usually a technical abbreviation) of the mesure, 
		as defined in the config file and recognize by this class. */
		string m_measureName;

		/** Pointer to a MeteoData's measure */
		double MeteoData::* m_measurePtr;

		};

		/**
		* @class FilterValue1Stn
		* @brief Filtering of a measure's values for single-station filters (abstract base class).
		* @author Florian Hof
		* @date   2009-03-19
		*/
		class FilterValue1Stn : public FilterBase1Stn {

	protected:

		// base definition

		FilterValue1Stn();

		virtual ~FilterValue1Stn();

		// check handling

		virtual void prepareCheck();

		// accessors

		/** 
		* Gets the name (usually a technical abbreviation) of the mesure, 
		* as defined in the config file and recognize by this class. 
		*/
		const string& getMeasureName() const;

		/**
		* Get the measure's value of a meteo data. 
		* The measure to return is determined by the corresponding filter's parameter.
		* @param data   The reference of the MeteoData from which the value will be taken.
		* @return The reference to the requested value; the value is thus editable. 
		*/
		double& getMeasureValue(MeteoData& data) const;
		const double& getMeasureValue(const MeteoData& data) const;

	private:

		// real implementation
		FilterValue* m_filterValue;

		};

		/**
		* @class FilterValueMStn
		* @brief Filtering of a measure's values for multi-stations filters (abstract base class).
		* @author Florian Hof
		* @date   2009-03-19
		*/
		class FilterValueMStn : public FilterBaseMStn {

	protected:

		// base definition

		FilterValueMStn();

		virtual ~FilterValueMStn();

		// check handling

		virtual void prepareCheck();

		// accessors

		/** 
		* Gets the name (usually a technical abbreviation) of the mesure, 
		* as defined in the config file and recognize by this class. 
		*/
		const string& getMeasureName() const;

		/**
		* Get the measure's value of a meteo data. 
		* The measure to return is determined by the corresponding filter's parameter.
		* @param data   The reference of the MeteoData from which the value will be taken.
		* @return The reference to the requested value; the value is thus editable. 
		*/
		double& getMeasureValue(MeteoData& data) const;
		const double& getMeasureValue(const MeteoData& data) const;

	private:

		// real implementation
		FilterValue* m_filterValue;

};

#endif 
