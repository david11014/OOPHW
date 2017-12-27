#include <string>
using namespace std;
#ifndef _TICK_TAC_TOE_BOARD_H
#define _TICK_TAC_TOE_BOARD_H
#define BOARD_SIZE 3
class TickTacToeBoard
{
public:
	
	int Board[BOARD_SIZE][BOARD_SIZE]; //盤面
	int NowPlayer; //現在玩家編號
	string PlayerName[2]; //玩家編號對應名字

	TickTacToeBoard();
	TickTacToeBoard(string P1, string P2);

	int SetBoard(int x, int y);
	int GetBoard(int x, int y);
	void SetPlayer(string P1, string P2);
	void Reset();
private:

	int CheckWiner();

};

#endif // !_TICK_TAC_TOE_BOARD_H





