#ifndef ViewH
#define ViewH

using namespace std;

#include <SFML/Graphics.hpp>
#include "RScore/Landscape.h"
#include "RScore/Community.h"

class BatchView {
public:
	BatchView(sf::RenderWindow& window, Landscape* pLand, Community* pCommunity);

	void collectUserInput(sf::RenderWindow& window);
	void drawLandscape(sf::RenderWindow& window);
	void drawCommunity(sf::RenderWindow& window, Species* pSpecies, const int& yr, const int& gen);

private:
	Landscape* pLandscape;
	Community* pComm;

	unsigned int cellSize;
	int dimX, dimY;

	const float relSizeLegend = 0.05; // height of time+gen label relative to dimX

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

	bool mustPause = false;
};

#endif ViewH