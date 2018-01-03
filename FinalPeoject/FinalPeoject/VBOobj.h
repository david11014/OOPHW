#pragma once

#ifndef __VBOobj_H__
#define __VBOobj_H__

class VBOobj
{
public:
	unsigned int* VBO;
	size_t VBOsize;

	VBOobj();
	VBOobj(const VBOobj& V);
	virtual ~VBOobj();

	VBOobj(size_t s);

	bool releaseVBO();
};

void genVBO(VBOobj* obj);
void deleteVBO(VBOobj* obj);

#endif