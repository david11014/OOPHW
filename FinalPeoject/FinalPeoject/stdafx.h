// stdafx.h : �i�b�����Y�ɤ��]�t�зǪ��t�� Include �ɡA
// �άO�g�`�ϥΫo�ܤ��ܧ�
// �M�ױM�� Include �ɮ�
#pragma once

// TODO: �b���Ѧұz���{���һݭn����L���Y

#define NOMINMAX

#ifdef _DEBUG

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#ifndef DBG_NEW
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
#define new DBG_NEW
#endif
#endif  // _DEBUG

#define _USE_MATH_DEFINES
#include "time_measure.h"

static mytime::mytime mtime;
static mytime::mytime mtime2;
static mytime::mytime mtime3;

#include <fstream>
#include <sstream>
