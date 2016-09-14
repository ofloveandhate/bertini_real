#include "data_type.hpp"
#include "programConfiguration.hpp"













int compare_integers_decreasing(const void * left_in, const void * right_in)
{
	
	int left = *(const int*)left_in;
	int right = *(const int*)right_in;
	
	
	if (left<right) {
		return 1;
	}
	else if (left>right){
		return -1;
	}
	else{
		return 0;
	}
	
}

int compare_integers_increasing(const void * left_in, const void * right_in)
{
	
	int left = *(const int*)left_in;
	int right = *(const int*)right_in;
	
	
	if (left>right) {
		return 1;
	}
	else if(left < right){
		return -1;
	}
	else{
		return 0;
	}
	
}













