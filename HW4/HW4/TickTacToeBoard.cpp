#include "TickTacToeBoard.h"



TickTacToeBoard::TickTacToeBoard()
{
	PlayerName[0] = "玩家1";
	PlayerName[1] = "玩家2";
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
			Board[i][j] = -1;
		}
	}
}

int TickTacToeBoard::CheckWiner()
{
	int winner = -1;

	// 斜下 
	winner = Board[0][0];
	for (int i = 0; i < BOARD_SIZE; i++)
	{	
		if (Board[0][0] != Board[i][i])
		{
			winner = -1;
			break;
		}
	}
	if (winner != -1)
		return winner;
	
	// 斜上 
	winner = Board[BOARD_SIZE - 1][0];
	for (int i = 0; i < BOARD_SIZE; i++)
	{
		if (Board[BOARD_SIZE - 1][0] != Board[BOARD_SIZE - i - 1][i])
		{
			winner = -1;
			break;
		}
	}
	if (winner != -1)
		return winner;

	//垂直
	for (int i = 0; i < BOARD_SIZE; i++)
	{
		winner = Board[i][0];

		for (int j = 0; j < BOARD_SIZE; j++)
		{
			if (Board[i][0] != Board[i][j])
			{
				winner = -1;
				break;
			}
		}

		if (winner != -1)
			return winner;
	}

	//水平
	for (int i = 0; i < BOARD_SIZE; i++)
	{
		winner = Board[0][i];

		for (int j = 0; j < BOARD_SIZE; j++)
		{
			if (Board[0][i] != Board[j][i])
			{
				winner = -1;
				break;
			}
		}

		if (winner != -1)
			return winner;
	}		

	//平手
	winner = -2;
	for (int i = 0; i < BOARD_SIZE; i++)
	{	
		for (int j = 0; j < BOARD_SIZE; j++)
		{
			if (Board[i][j] == -1)
			{
				winner = -1;
			}
		}
	}

	return winner;
}
