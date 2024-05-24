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
	void drawCommunity(sf::RenderWindow& window, Species* pSpecies);

private:
	Landscape* pLandscape;
	Community* pComm;

	unsigned int cellSize;
	int dimX, dimY;

	float indRadius = 1.0;
	sf::CircleShape indShape = sf::CircleShape(indRadius);
	sf::Color indColour = sf::Color::Red;

	const vector <sf::Color> habitatPalette{
	sf::Color::White,
	sf::Color::Blue,
	sf::Color::Cyan,
	sf::Color::Green,
	sf::Color::Yellow,
	sf::Color::Black,
	sf::Color::Magenta
	};
};

#endif ViewH