#pragma once
#include<regex>
#include<string>
#include "TickTacToeBoard.h"
namespace HW4 {

	using namespace std;
	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// HW4Form 的摘要
	/// </summary>
	public ref class HW4Form : public System::Windows::Forms::Form
	{
	public:
		TickTacToeBoard *GameBoard;

		HW4Form(void)
		{
			InitializeComponent();
			//
			//TODO:  在此加入建構函式程式碼
			//
			GameBoard = new TickTacToeBoard();
		}

	protected:
		/// <summary>
		/// 清除任何使用中的資源。
		/// </summary>
		~HW4Form()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel1;
	protected:
	private: System::Windows::Forms::StatusStrip^  statusStrip1;
	private: System::Windows::Forms::ToolStripStatusLabel^  NowPlayerLabel;

	private: System::Windows::Forms::GroupBox^  groupBox1;
	private: System::Windows::Forms::RadioButton^  P2RadioButton;
	private: System::Windows::Forms::RadioButton^  P1RadioButton;
	private: System::Windows::Forms::PictureBox^  pictureBox3_3;
	private: System::Windows::Forms::PictureBox^  pictureBox3_2;
	private: System::Windows::Forms::PictureBox^  pictureBox3_1;
	private: System::Windows::Forms::PictureBox^  pictureBox2_3;
	private: System::Windows::Forms::PictureBox^  pictureBox2_1;
	private: System::Windows::Forms::PictureBox^  pictureBox1_3;
	private: System::Windows::Forms::PictureBox^  pictureBox1_2;
	private: System::Windows::Forms::PictureBox^  pictureBox1_1;
	private: System::Windows::Forms::PictureBox^  pictureBox2_2;
	private: System::Windows::Forms::Button^  ControlButton;


	private: System::Windows::Forms::Label^  TimeLabel;
	private: System::Windows::Forms::ToolStripStatusLabel^  MouseLocLecLabel;
	private: System::Windows::Forms::Timer^  timer1;
	private: System::ComponentModel::IContainer^  components;

	private:
		/// <summary>
		/// 設計工具所需的變數。
		/// </summary>


#pragma region Windows Form Designer generated code
		/// <summary>
		/// 此為設計工具支援所需的方法 - 請勿使用程式碼編輯器修改
		/// 這個方法的內容。
		/// </summary>
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			this->tableLayoutPanel1 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->pictureBox3_3 = (gcnew System::Windows::Forms::PictureBox());
			this->pictureBox3_2 = (gcnew System::Windows::Forms::PictureBox());
			this->pictureBox3_1 = (gcnew System::Windows::Forms::PictureBox());
			this->pictureBox1_2 = (gcnew System::Windows::Forms::PictureBox());
			this->pictureBox2_2 = (gcnew System::Windows::Forms::PictureBox());
			this->pictureBox1_3 = (gcnew System::Windows::Forms::PictureBox());
			this->pictureBox2_3 = (gcnew System::Windows::Forms::PictureBox());
			this->pictureBox2_1 = (gcnew System::Windows::Forms::PictureBox());
			this->pictureBox1_1 = (gcnew System::Windows::Forms::PictureBox());
			this->statusStrip1 = (gcnew System::Windows::Forms::StatusStrip());
			this->NowPlayerLabel = (gcnew System::Windows::Forms::ToolStripStatusLabel());
			this->MouseLocLecLabel = (gcnew System::Windows::Forms::ToolStripStatusLabel());
			this->groupBox1 = (gcnew System::Windows::Forms::GroupBox());
			this->P2RadioButton = (gcnew System::Windows::Forms::RadioButton());
			this->P1RadioButton = (gcnew System::Windows::Forms::RadioButton());
			this->ControlButton = (gcnew System::Windows::Forms::Button());
			this->TimeLabel = (gcnew System::Windows::Forms::Label());
			this->timer1 = (gcnew System::Windows::Forms::Timer(this->components));
			this->tableLayoutPanel1->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox3_3))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox3_2))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox3_1))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1_2))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox2_2))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1_3))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox2_3))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox2_1))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1_1))->BeginInit();
			this->statusStrip1->SuspendLayout();
			this->groupBox1->SuspendLayout();
			this->SuspendLayout();
			// 
			// tableLayoutPanel1
			// 
			this->tableLayoutPanel1->CellBorderStyle = System::Windows::Forms::TableLayoutPanelCellBorderStyle::Single;
			this->tableLayoutPanel1->ColumnCount = 3;
			this->tableLayoutPanel1->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent,
				33.33333F)));
			this->tableLayoutPanel1->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent,
				33.33333F)));
			this->tableLayoutPanel1->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent,
				33.33333F)));
			this->tableLayoutPanel1->Controls->Add(this->pictureBox3_3, 2, 2);
			this->tableLayoutPanel1->Controls->Add(this->pictureBox3_2, 1, 2);
			this->tableLayoutPanel1->Controls->Add(this->pictureBox3_1, 0, 2);
			this->tableLayoutPanel1->Controls->Add(this->pictureBox1_2, 1, 0);
			this->tableLayoutPanel1->Controls->Add(this->pictureBox2_2, 1, 1);
			this->tableLayoutPanel1->Controls->Add(this->pictureBox2_1, 0, 1);
			this->tableLayoutPanel1->Controls->Add(this->pictureBox1_1, 0, 0);
			this->tableLayoutPanel1->Controls->Add(this->pictureBox2_3, 2, 1);
			this->tableLayoutPanel1->Controls->Add(this->pictureBox1_3, 2, 0);
			this->tableLayoutPanel1->Location = System::Drawing::Point(12, 12);
			this->tableLayoutPanel1->Margin = System::Windows::Forms::Padding(0);
			this->tableLayoutPanel1->Name = L"tableLayoutPanel1";
			this->tableLayoutPanel1->RowCount = 3;
			this->tableLayoutPanel1->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 33.33333F)));
			this->tableLayoutPanel1->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 33.33333F)));
			this->tableLayoutPanel1->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 33.33333F)));
			this->tableLayoutPanel1->Size = System::Drawing::Size(478, 488);
			this->tableLayoutPanel1->TabIndex = 0;
			// 
			// pictureBox3_3
			// 
			this->pictureBox3_3->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->pictureBox3_3->Location = System::Drawing::Point(319, 325);
			this->pictureBox3_3->Margin = System::Windows::Forms::Padding(0);
			this->pictureBox3_3->Name = L"pictureBox3_3";
			this->pictureBox3_3->Size = System::Drawing::Size(158, 162);
			this->pictureBox3_3->TabIndex = 8;
			this->pictureBox3_3->TabStop = false;
			this->pictureBox3_3->Tag = L"33";
			this->pictureBox3_3->Click += gcnew System::EventHandler(this, &HW4Form::Board_Click);
			this->pictureBox3_3->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &HW4Form::Board_MouseMove);
			// 
			// pictureBox3_2
			// 
			this->pictureBox3_2->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->pictureBox3_2->Location = System::Drawing::Point(160, 325);
			this->pictureBox3_2->Margin = System::Windows::Forms::Padding(0);
			this->pictureBox3_2->Name = L"pictureBox3_2";
			this->pictureBox3_2->Size = System::Drawing::Size(158, 162);
			this->pictureBox3_2->TabIndex = 7;
			this->pictureBox3_2->TabStop = false;
			this->pictureBox3_2->Tag = L"32";
			this->pictureBox3_2->Click += gcnew System::EventHandler(this, &HW4Form::Board_Click);
			this->pictureBox3_2->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &HW4Form::Board_MouseMove);
			// 
			// pictureBox3_1
			// 
			this->pictureBox3_1->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->pictureBox3_1->Location = System::Drawing::Point(1, 325);
			this->pictureBox3_1->Margin = System::Windows::Forms::Padding(0);
			this->pictureBox3_1->Name = L"pictureBox3_1";
			this->pictureBox3_1->Size = System::Drawing::Size(158, 162);
			this->pictureBox3_1->TabIndex = 6;
			this->pictureBox3_1->TabStop = false;
			this->pictureBox3_1->Tag = L"31";
			this->pictureBox3_1->Click += gcnew System::EventHandler(this, &HW4Form::Board_Click);
			this->pictureBox3_1->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &HW4Form::Board_MouseMove);
			// 
			// pictureBox1_2
			// 
			this->pictureBox1_2->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->pictureBox1_2->Location = System::Drawing::Point(160, 1);
			this->pictureBox1_2->Margin = System::Windows::Forms::Padding(0);
			this->pictureBox1_2->Name = L"pictureBox1_2";
			this->pictureBox1_2->Size = System::Drawing::Size(158, 161);
			this->pictureBox1_2->TabIndex = 1;
			this->pictureBox1_2->TabStop = false;
			this->pictureBox1_2->Tag = L"12";
			this->pictureBox1_2->Click += gcnew System::EventHandler(this, &HW4Form::Board_Click);
			this->pictureBox1_2->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &HW4Form::Board_MouseMove);
			// 
			// pictureBox2_2
			// 
			this->pictureBox2_2->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->pictureBox2_2->Location = System::Drawing::Point(160, 163);
			this->pictureBox2_2->Margin = System::Windows::Forms::Padding(0);
			this->pictureBox2_2->Name = L"pictureBox2_2";
			this->pictureBox2_2->Size = System::Drawing::Size(158, 161);
			this->pictureBox2_2->TabIndex = 4;
			this->pictureBox2_2->TabStop = false;
			this->pictureBox2_2->Tag = L"22";
			this->pictureBox2_2->Click += gcnew System::EventHandler(this, &HW4Form::Board_Click);
			this->pictureBox2_2->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &HW4Form::Board_MouseMove);
			// 
			// pictureBox1_3
			// 
			this->pictureBox1_3->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->pictureBox1_3->Location = System::Drawing::Point(319, 1);
			this->pictureBox1_3->Margin = System::Windows::Forms::Padding(0);
			this->pictureBox1_3->Name = L"pictureBox1_3";
			this->pictureBox1_3->Size = System::Drawing::Size(158, 161);
			this->pictureBox1_3->TabIndex = 2;
			this->pictureBox1_3->TabStop = false;
			this->pictureBox1_3->Tag = L"13";
			this->pictureBox1_3->Click += gcnew System::EventHandler(this, &HW4Form::Board_Click);
			this->pictureBox1_3->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &HW4Form::Board_MouseMove);
			// 
			// pictureBox2_3
			// 
			this->pictureBox2_3->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->pictureBox2_3->Location = System::Drawing::Point(319, 163);
			this->pictureBox2_3->Margin = System::Windows::Forms::Padding(0);
			this->pictureBox2_3->Name = L"pictureBox2_3";
			this->pictureBox2_3->Size = System::Drawing::Size(158, 161);
			this->pictureBox2_3->TabIndex = 5;
			this->pictureBox2_3->TabStop = false;
			this->pictureBox2_3->Tag = L"23";
			this->pictureBox2_3->Click += gcnew System::EventHandler(this, &HW4Form::Board_Click);
			this->pictureBox2_3->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &HW4Form::Board_MouseMove);
			// 
			// pictureBox2_1
			// 
			this->pictureBox2_1->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->pictureBox2_1->Location = System::Drawing::Point(1, 163);
			this->pictureBox2_1->Margin = System::Windows::Forms::Padding(0);
			this->pictureBox2_1->Name = L"pictureBox2_1";
			this->pictureBox2_1->Size = System::Drawing::Size(158, 161);
			this->pictureBox2_1->TabIndex = 3;
			this->pictureBox2_1->TabStop = false;
			this->pictureBox2_1->Tag = L"21";
			this->pictureBox2_1->Click += gcnew System::EventHandler(this, &HW4Form::Board_Click);
			this->pictureBox2_1->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &HW4Form::Board_MouseMove);
			// 
			// pictureBox1_1
			// 
			this->pictureBox1_1->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->pictureBox1_1->Location = System::Drawing::Point(1, 1);
			this->pictureBox1_1->Margin = System::Windows::Forms::Padding(0);
			this->pictureBox1_1->Name = L"pictureBox1_1";
			this->pictureBox1_1->Size = System::Drawing::Size(158, 161);
			this->pictureBox1_1->TabIndex = 0;
			this->pictureBox1_1->TabStop = false;
			this->pictureBox1_1->Tag = L"11";
			this->pictureBox1_1->Click += gcnew System::EventHandler(this, &HW4Form::Board_Click);
			this->pictureBox1_1->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &HW4Form::Board_MouseMove);
			// 
			// statusStrip1
			// 
			this->statusStrip1->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(2) {
				this->NowPlayerLabel,
					this->MouseLocLecLabel
			});
			this->statusStrip1->Location = System::Drawing::Point(0, 519);
			this->statusStrip1->Name = L"statusStrip1";
			this->statusStrip1->Size = System::Drawing::Size(779, 22);
			this->statusStrip1->TabIndex = 1;
			this->statusStrip1->Text = L"statusStrip1";
			// 
			// NowPlayerLabel
			// 
			this->NowPlayerLabel->Name = L"NowPlayerLabel";
			this->NowPlayerLabel->Size = System::Drawing::Size(38, 17);
			this->NowPlayerLabel->Text = L"玩家1";
			// 
			// MouseLocLecLabel
			// 
			this->MouseLocLecLabel->Name = L"MouseLocLecLabel";
			this->MouseLocLecLabel->Size = System::Drawing::Size(90, 17);
			this->MouseLocLecLabel->Text = L"滑鼠目前位置(,)";
			// 
			// groupBox1
			// 
			this->groupBox1->Controls->Add(this->P2RadioButton);
			this->groupBox1->Controls->Add(this->P1RadioButton);
			this->groupBox1->Font = (gcnew System::Drawing::Font(L"新細明體", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(136)));
			this->groupBox1->Location = System::Drawing::Point(500, 12);
			this->groupBox1->Name = L"groupBox1";
			this->groupBox1->Size = System::Drawing::Size(256, 90);
			this->groupBox1->TabIndex = 2;
			this->groupBox1->TabStop = false;
			this->groupBox1->Text = L"先手(O)";
			// 
			// P2RadioButton
			// 
			this->P2RadioButton->AutoSize = true;
			this->P2RadioButton->Font = (gcnew System::Drawing::Font(L"新細明體", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(136)));
			this->P2RadioButton->Location = System::Drawing::Point(7, 51);
			this->P2RadioButton->Name = L"P2RadioButton";
			this->P2RadioButton->Size = System::Drawing::Size(74, 23);
			this->P2RadioButton->TabIndex = 1;
			this->P2RadioButton->Text = L"玩家2";
			this->P2RadioButton->UseVisualStyleBackColor = true;
			// 
			// P1RadioButton
			// 
			this->P1RadioButton->AutoSize = true;
			this->P1RadioButton->Checked = true;
			this->P1RadioButton->Font = (gcnew System::Drawing::Font(L"新細明體", 14.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(136)));
			this->P1RadioButton->Location = System::Drawing::Point(7, 22);
			this->P1RadioButton->Name = L"P1RadioButton";
			this->P1RadioButton->Size = System::Drawing::Size(74, 23);
			this->P1RadioButton->TabIndex = 0;
			this->P1RadioButton->TabStop = true;
			this->P1RadioButton->Text = L"玩家1";
			this->P1RadioButton->UseVisualStyleBackColor = true;
			// 
			// ControlButton
			// 
			this->ControlButton->Font = (gcnew System::Drawing::Font(L"新細明體", 18, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(136)));
			this->ControlButton->Location = System::Drawing::Point(500, 382);
			this->ControlButton->Name = L"ControlButton";
			this->ControlButton->Size = System::Drawing::Size(251, 67);
			this->ControlButton->TabIndex = 3;
			this->ControlButton->Text = L"Start";
			this->ControlButton->UseVisualStyleBackColor = true;
			this->ControlButton->Click += gcnew System::EventHandler(this, &HW4Form::ControlButton_Click);
			// 
			// TimeLabel
			// 
			this->TimeLabel->AutoSize = true;
			this->TimeLabel->Font = (gcnew System::Drawing::Font(L"新細明體", 27.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(136)));
			this->TimeLabel->Location = System::Drawing::Point(500, 463);
			this->TimeLabel->Name = L"TimeLabel";
			this->TimeLabel->Size = System::Drawing::Size(182, 37);
			this->TimeLabel->TabIndex = 4;
			this->TimeLabel->Text = L"00:00:00.00";
			// 
			// timer1
			// 
			this->timer1->Interval = 10;
			// 
			// HW4Form
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 12);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(779, 541);
			this->Controls->Add(this->TimeLabel);
			this->Controls->Add(this->ControlButton);
			this->Controls->Add(this->groupBox1);
			this->Controls->Add(this->statusStrip1);
			this->Controls->Add(this->tableLayoutPanel1);
			this->Name = L"HW4Form";
			this->Text = L"HW4Form";
			this->tableLayoutPanel1->ResumeLayout(false);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox3_3))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox3_2))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox3_1))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1_2))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox2_2))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1_3))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox2_3))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox2_1))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1_1))->EndInit();
			this->statusStrip1->ResumeLayout(false);
			this->statusStrip1->PerformLayout();
			this->groupBox1->ResumeLayout(false);
			this->groupBox1->PerformLayout();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion

private:
	void ResetPictureBoxColor()
	{		
		pictureBox3_3->CreateGraphics()->Clear(Control::DefaultBackColor);
		pictureBox3_2->CreateGraphics()->Clear(Control::DefaultBackColor);
		pictureBox3_1->CreateGraphics()->Clear(Control::DefaultBackColor);
		pictureBox2_3->CreateGraphics()->Clear(Control::DefaultBackColor);
		pictureBox2_1->CreateGraphics()->Clear(Control::DefaultBackColor);
		pictureBox1_3->CreateGraphics()->Clear(Control::DefaultBackColor);
		pictureBox1_2->CreateGraphics()->Clear(Control::DefaultBackColor);
		pictureBox1_1->CreateGraphics()->Clear(Control::DefaultBackColor);
		pictureBox2_2->CreateGraphics()->Clear(Control::DefaultBackColor);
	}
	void StartGame()
	{
		ControlButton->Text = "Stop";
		GameBoard->NowPlayer = 0;
		if (P1RadioButton->Checked == true)
		{
			GameBoard->SetPlayer("玩家1", "玩家2");
			P2RadioButton->Enabled = false;
		}
		if (P2RadioButton->Checked == true)
		{
			GameBoard->SetPlayer("玩家2", "玩家1");
			P1RadioButton->Enabled = false;
		}		
		
		NowPlayerLabel->Text = gcnew String((GameBoard->PlayerName[GameBoard->NowPlayer]).c_str());
		timer1->Start();
	}
	void EndGame()
	{
		ResetPictureBoxColor();
		ControlButton->Text = "Start";

		P1RadioButton->Enabled = true;
		P2RadioButton->Enabled = true;
		timer1->Stop();
	}

private: System::Void Board_MouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) {
	
	System::Windows::Forms::PictureBox^ nowPB = (System::Windows::Forms::PictureBox^)sender;
	
	int tag_num = Convert::ToInt32(nowPB->Tag);
	int x = tag_num / 10;
	int y = tag_num % 10;	

	MouseLocLecLabel->Text = "滑鼠目前位置: (" + x.ToString() + ", " + y.ToString() + ")";
	
}
private: System::Void Board_Click(System::Object^  sender, System::EventArgs^  e) {

	if (timer1->Enabled == true)
	{
		System::Windows::Forms::PictureBox^ nowPB = (System::Windows::Forms::PictureBox^)sender;

		Graphics ^G = nowPB->CreateGraphics();
		float margin = nowPB->Width*0.1;
		Pen penX(Color::Black, 5);
		Pen penO(Color::Red, 5);

		if (GameBoard->NowPlayer == 0)
		{
			G->DrawEllipse(%penO, margin, margin, nowPB->Width - margin * 2, nowPB->Height - margin * 2);
		}
		else
		{
			G->DrawLine(%penX, margin, margin, nowPB->Width - margin, nowPB->Height - margin);
			G->DrawLine(%penX, nowPB->Width - margin, margin, margin, nowPB->Height - margin);
		}

		int tag_num = Convert::ToInt32(nowPB->Tag);
		int x = tag_num / 10;
		int y = tag_num % 10;

		GameBoard->SetBoard(x - 1, y - 1);		
		NowPlayerLabel->Text = gcnew String((GameBoard->PlayerName[GameBoard->NowPlayer]).c_str());
	}
	else
		return;

}
private: System::Void ControlButton_Click(System::Object^  sender, System::EventArgs^  e) {

	if (timer1->Enabled == false)
		StartGame();
	else
		EndGame();		
}
};
}
