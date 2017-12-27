#include "TickTacToeBoard.h"



TickTacToeBoard::TickTacToeBoard()
{
	PlayerName[0] = "ª±®a1";
	PlayerName[1] = "ª±®a2";
	NowPlayer = 0;
	Reset();
}

TickTacToeBoard::TickTacToeBoard(string P1,string P2)
{
	SetPlayer(P1, P2);
	NowPlayer = 0;
	Reset();
}

int TickTacToeBoard::SetBoard(int x, int y)
{
	Board[x][y] = NowPlayer;
	NowPlayer = !NowPlayer;
	return CheckWiner();
}

int TickTacToeBoard::GetBoard(int x, int y)
{
	return Board[x][y];
}

void TickTacToeBoard::SetPlayer(string P1, string P2)
{
	PlayerName[0] = P1;
	PlayerName[1] = P2;
}

void TickTacToeBoard::Reset()
{
	NowPlayer = 0;
	for (int i = 0; i < BOARD_SIZE; i++)
	{
		for (int j= 0; j < BOARD_SIZE; j++)
		{
			Board[i][j] == -1;
		}
	}
}

int TickTacToeBoard::CheckWiner()
{
	return -1;
}
