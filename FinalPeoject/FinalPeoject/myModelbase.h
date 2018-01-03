#pragma once

#ifndef __myModelbase_H__
#define __myModelbase_H__

template<class T>
using sp = std::shared_ptr<T>;
template<class T>
using up = std::unique_ptr<T>;
template<class T>
using vectors = std::vector<std::vector<T>>;
template<class T>
using spvector = std::vector<sp<T>>;
template<class T>
using spvectors = std::vector<spvector<T>>;

namespace myModel
{
	typedef size_t model_ind_t;
	typedef size_t tri_ind_t;
	typedef unsigned char vertex_ind_t;
	typedef vertex_ind_t edge_ind_t;

	typedef vector<tri_ind_t> inddata;
	typedef vector<inddata> inddatas;

	typedef unsigned int gl_tri_ind_t;
}

#endif