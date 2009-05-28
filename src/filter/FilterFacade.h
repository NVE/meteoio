#ifndef FILTERFACADE_H_INCLUDED
#define FILTERFACADE_H_INCLUDED

#include "FilterBase1Stn.h"
#include "FilterBaseMStn.h"

using namespace std;

/**
 * @class FilterFacade
 * @brief Facade for the handling and running of filters.
 * @author Florian Hof
 * @date   2009-03-02
 */
class FilterFacade {

	public:

		// base definition

		/** 
		* Default constructor. 
		* Does only construct the object; defining filters and all the stuff has to be done manually. 
		*/
		FilterFacade();

		/** 
		* Quick and easy constructor. 
		* Its call replaces calls to FilterFacade(), addActiveFilters(filename), prepareCheck(), 
		* and also clearActiveFilters() when destructing. 
		* @param filename   [in] The file path of the .ini config file that defines the parameters. 
		*/
		FilterFacade(const string& filename);

		/** 
		* Destructor. 
		* Note that the filters are not deleted here, unless you use the quick and easy constructor. 
		* If filters are not needed anymore, call clearActiveFilters before. 
		*/
		~FilterFacade();

		// types of filters handling

		/** 
		* Add a type of single-station filter to the known filters directory. 
		* Every filter has to be registered for being available for use. The registration as to be done
		* 1) either in the implementation of FilterFacade::registerFilters();
		* 2) or elsewhere where its get executed, such as in the main function; 
		*    note that the filter's implementation is not such a place. 
		* @param name   The unique name of the filter's type (usually its class name). 
		* @param creator   A function that returns a new instance of the filter's type. 
		*/
		static void registerFilter(const string& name, const Filter1StnCreator creator);

		/** 
		* Add a type of multi-stations filter to the known filters directory. 
		* Every filter has to be registered for being available for use. The registration as to be done
		* 1) either in the implementation of FilterFacade::registerFilters();
		* 2) or elsewhere where its get executed, such as in the main function; 
		*    note that the filter's implementation is not such a place. 
		* @param name   The unique name of the filter's type (usually its class name). 
		* @param creator   A function that returns a new instance of the filter's type. 
		*/
		static void registerFilter(const string& name, const FilterMStnCreator creator);

		/**
		* Fill the list of known types of filter.
		* This is a convenient place for listing all known filters. 
		*/
		static void registerAllFilters();

		/** 
		* Lists the known types of single-station filters. 
		* @return The directory of known types of filters, mapping each filter's name to the corresponding filter's creator. 
		*/
		static map<string, Filter1StnCreator>* getRegisteredFilters1Stn();

		/** 
		* Lists the known types of multi-stations filters. 
		* @return The directory of known types of filters, mapping each filter's name to the corresponding filter's creator. 
		*/
		static map<string, FilterMStnCreator>* getRegisteredFiltersMStn();

		// filters handling

		/**
		* Get the active single-station filters (i.e. the check-list)
		*/
		vector<FilterBase1Stn*>* getActiveFilters1Stn();

		/**
		* Get the active multi-stations filters (i.e. the check-list)
		*/
		vector<FilterBaseMStn*>* getActiveFiltersMStn();

		/**
		* Add filters to the check-lists, readen from a .INI config file. 
		* @param filename   [in] The file path of the .ini config file that defines the filters and their parameters. 
		*/
		void addActiveFilters(const string& filename);

		/** 
		* Add a single-station filter to the check-list.
		* @param filter   The reference to the filter to add. 
		*/
		void addActiveFilter(FilterBase1Stn* filter);

		/** 
		* Add a multi-stations filter to the check-list.
		* @param filter   The reference to the filter to add. 
		*/
		void addActiveFilter(FilterBaseMStn* filter);

		/** 
		* Remove a single-station filter from the check-list.
		* Note that the filter is not deleted here; if not needed anymore, the filter has to be freeed and its memory released. 
		* @param filter   The reference to the filter to remove. 
		* @return True if the filter was successfully deleted, false otherwise (not found, etc.). 
		*/ 
		bool removeActiveFilter(const FilterBase1Stn* filter);

		/** 
		* Remove a multi-stations filter from the check-list.
		* Note that the filter is not deleted here; if not needed anymore, the filter has to be freeed and its memory released. 
		* @param filter   The reference to the filter to remove. 
		* @return True if the filter was successfully deleted, false otherwise (not found, etc.). 
		*/ 
		bool removeActiveFilter(const FilterBaseMStn* filter);

		/** 
		* Remove all single-station and multi-stations filters from the check-lists.
		* Note that the filters are not deleted here; if not needed anymore, call clearActiveFilters instead. 
		*/ 
		void removeActiveFilters();

		/** 
		* Remove and delete all single-station and multi-stations filters from the check-lists.
		* Note that the filters are deleted here and cannot be used anymore; its memory released. 
		*/ 
		void clearActiveFilters();

		// check handling

		/**
		* Returns the minimal size of the window that the filters requires.
		* Both conditions has to be fullfilled for the filter to make the proper decision (otherwise, a warning is printed during doCheck). 
		* Before calling getMinimalWindow, the parameters have to be filled and the prepareCheck method called!
		* @param minNbPoints   [out] The minimal number of measures. Returns 1 for single-value filters (for ex. min-max).
		* @param minDeltaTime   [out] The minimal delta of the time frame. Returns an empty delta for single-value filters. 
		*/
		void getMinimalWindow(unsigned int& minNbPoints, Date_IO& minDeltaTime);

		/**
		* Make some initialisation and parameters' validation before the real checks.
		* Has to be called before any doCheck and getMinimalWindow methods and after any filters' modification
		*/
		void prepareCheck();

		/**
		* Check the meteo data at a single time with every active single-station filters.
		* Before calling doCheck, the filter's parameters have to be filled and the prepareCheck method called. 
		* @param unfilteredMeteoBuffer   [in] The buffer of unfiltered meteo data.
		* @param filteredMeteoBuffer   [in out] The buffer of filtered meteo data. May content already-filtered 
		*                              data at some previous time or at the time to check, or may also be empty. 
		* @param iUnfilteredElement   [in] The index in the unfilteredMeteoBuffer of the meteo data to filter, 
		*                             thus specifying the date and time. 
		*/
		void doCheck(MeteoBuffer& unfilteredMeteoBuffer, MeteoBuffer& filteredMeteoBuffer, unsigned int iUnfilteredElement);

		/**
		* Check the meteo data at a single time with every active single-station and multi-stations filters.
		* Typically used on ressampled meteodata of several stations. 
		* Before calling doCheck, the filter's parameters have to be filled and the prepareCheck method called. 
		* @param stations   [in] The list of stations to which the meteo data corresponds.
		* @param meteoBuffers   [in] For each station, the buffer (over the past) of meteo data, used by filters that needs a time frame. 
		* @param meteoDatas   [in out] For each station, the meteo data to filter. 
		*/
		void doCheck(vector<StationData>& stations, vector<MeteoBuffer>& meteoBuffers, vector<MeteoData>& meteoDatas);

	protected:

		/** directory of all known types of single-station filters, with their name and creator */
		static map<string, Filter1StnCreator> s_registeredFilters1StnByName;

		/** directory of all known types of multi-stations filters, with their name and creator */
		static map<string, FilterMStnCreator> s_registeredFiltersMStnByName;

		/** list of all active single-station filters (i.e. check-list) */
		vector<FilterBase1Stn*> m_activeFilters1Stn;

		/** list of all active multi-stations filters (i.e. check-list) */
		vector<FilterBaseMStn*> m_activeFiltersMStn;

	private: 

		/** States if the quick and easy constructor was used, which implies to call clearActiveFilters at destruction */
		bool m_isEasyMode;

};
  
#endif
