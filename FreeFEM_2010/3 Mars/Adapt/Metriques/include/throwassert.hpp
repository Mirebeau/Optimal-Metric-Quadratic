// -*- Mode : c++ -*-
//
// SUMMARY  :      
// USAGE    :        
// ORG      : 
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

/*
 
 This file is part of Freefem++
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef THROWASSERT
#define THROWASSERT
#include <iostream>

//#ifdef __INTEL__
#define cerr cout 
//#endif

#include "error.hpp"

#ifdef NDEBUG
#define throwassert(i)  ( (void) 0)
#else
#define throwassert(condition)  ((condition) ? ( (void) 0) : throw(ErrorAssert(#condition,__FILE__, __LINE__)))
 
#undef assert
#define assert(condition) throwassert(condition)
#endif
// an unremovable assert : ffassert
#define ffassert(condition)  ((condition) ? ( (void) 0) : throw(ErrorAssert(#condition,__FILE__, __LINE__)))
#define InternalError(message) throw(ErrorInternal(message,__LINE__,__FILE__))
#endif
