// stdafx.h : 可在此標頭檔中包含標準的系統 Include 檔，
// 或是經常使用卻很少變更的
// 專案專用 Include 檔案
#pragma once

// TODO: 在此參考您的程式所需要的其他標頭

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
