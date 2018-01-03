#include "MainForm.h"

using namespace System;
using namespace System::Windows::Forms;
[STAThread]

void main()
{
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false);

	FinalPeoject::MainForm mainform;

	Application::Run(%mainform);

}