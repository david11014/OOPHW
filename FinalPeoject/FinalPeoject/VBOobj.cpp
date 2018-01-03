#include "stdafx.h"

//#include "gl/GLee.h"

#include "VBOobj.h"

VBOobj::VBOobj() :VBOsize(2), VBO(nullptr)
{
}

VBOobj::VBOobj(const VBOobj& V) :VBOsize(V.VBOsize), VBO(nullptr)
{
}

VBOobj::~VBOobj()
{
	if (VBO)
	{
		deleteVBO(this);
		delete[] VBO;
	}
}

VBOobj::VBOobj(size_t s) :VBOsize(s), VBO(nullptr)
{
}

bool VBOobj::releaseVBO()
{
	if (VBO)
	{
		deleteVBO(this);
		delete[] VBO;
		VBO = nullptr;
		return true;
	}
	return false;
}

//void genVBO(VBOobj* obj)
//{
//	if (!obj->VBO)
//	{
//		obj->VBO = new unsigned int[obj->VBOsize];
//		glGenBuffers((GLsizei)obj->VBOsize, obj->VBO);
//	}
//}
//
//void deleteVBO(VBOobj* obj)
//{
//	if (obj != nullptr)
//	{
//		if (glIsBuffer(obj->VBO[0]))
//		{
//			glDeleteBuffers((GLsizei)obj->VBOsize, obj->VBO);
//		}
//		else
//		{
//			std::cout << "release VBO fail.\n";
//		}
//	}
//}