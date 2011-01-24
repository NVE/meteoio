#include <meteoio/FilterBlock.h>

namespace mio {

FilterBlock::FilterBlock(const std::string& filter_name) : ProcessingBlock(filter_name) 
{

}

FilterBlock::~FilterBlock() {}

bool FilterBlock::is_soft(std::vector<std::string>& vec_args)
{
	if (vec_args.size() > 0){
		if (vec_args[0] == "soft"){
			vec_args.erase(vec_args.begin());
			return true;
		}
	}
	
	return false;
}

unsigned int FilterBlock::get_orientation(std::vector<std::string>& vec_args)
{
	if (vec_args.size() > 0){
		if (vec_args[0] == "left"){
			vec_args.erase(vec_args.begin());
			return FilterBlock::left;
		} else if (vec_args[0] == "right"){
			vec_args.erase(vec_args.begin());
			return FilterBlock::right;
		} else if (vec_args[0] == "center"){
			vec_args.erase(vec_args.begin());
			return FilterBlock::center;
		}
	}
	
	return FilterBlock::center; //the default
}

void FilterBlock::convert_args(const unsigned int& min_nargs, const unsigned int& max_nargs,
						 const std::vector<std::string>& vec_args, std::vector<double>& dbl_args)
{
	if ((vec_args.size() < min_nargs) || (vec_args.size() > max_nargs))
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName(), AT); 

	for (unsigned int ii=0; ii<vec_args.size(); ii++){
		double tmp = IOUtils::nodata;
		IOUtils::convertString(tmp, vec_args[ii]);
		dbl_args.push_back(tmp);
	}
}

}
