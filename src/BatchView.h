#ifndef ViewH
#define ViewH

using namespace std;

#include <SFML/Graphics.hpp>
#include "RScore/Landscape.h"
#include "RScore/Community.h"

class BatchView {
public:
	BatchView::BatchView(sf::RenderWindow& window, Landscape* pLand, Community* pCommunity);

	void collectUserInput(sf::RenderWindow& window);
	void drawLandscape(sf::RenderWindow& window);
	void drawCommunity(sf::RenderWindow& window);

private:
	Landscape* pLandscape;
	Community* pComm;

	unsigned int cellSize;
	int dimX, dimY;

	sf::Color cellColour = sf::Color(85, 51, 255); // Han purple

	const vector <sf::Color> habitatPalette{
	sf::Color::White,
	sf::Color::Blue,
	sf::Color::Red,
	sf::Color::Green,
	sf::Color::Yellow,
	sf::Color::Black,
	sf::Color::Magenta,
	sf::Color::Cyan
	};
};

#endif ViewH