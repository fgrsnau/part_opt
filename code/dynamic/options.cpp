#include "options.h"

template<> toption<bool>::operator bool()const{
	return *val!=0;
};