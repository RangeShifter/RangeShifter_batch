#ifndef ViewH
#define ViewH

using namespace std;

#include <SFML/Graphics.hpp>
#include "RScore/Landscape.h"
#include "RScore/Community.h"
#include <filesystem>

constexpr unsigned int dfltWinWidth = 1920u;
constexpr unsigned int dfltWinHeight = 1080u;
const float relSizeLegend = 0.10; // height of time+gen label relative to dimX

class BatchView {
public:
	BatchView(sf::RenderWindow& window, Landscape* pLand, Community* pCommunity,
		const int& maxYear, const int& maxGen);

	void collectUserInput(sf::RenderWindow& window);
	bool isPaused() const { return paused; }
	void pause() { paused = true; }
	void runPauseLoop(sf::RenderWindow& window);

	void drawLandscape(sf::RenderWindow& window);
	void drawCommunity(sf::RenderWindow& window, Species* pSpecies, const int& yr, const int& gen);

private:
	Landscape* pLandscape;
	Community* pComm;

	int maxYear, maxSeason;

	unsigned int cellSize;
	unsigned int winWidth, winHeight;

	sf::Font font;

	sf::Text timeLegendTxt, txtPaused;
	sf::RectangleShape timeLegendBg, txtPausedBg;

	float indRadius = 1.0;
	sf::CircleShape indShape = sf::CircleShape(indRadius);
	sf::Color indColour = sf::Color(255, 0, 110);

	const vector <sf::Color> habitatPalette{
	sf::Color(33, 158, 188),
	sf::Color(251, 133, 0),
	sf::Color(142, 202, 230),
	sf::Color::Green,
	sf::Color(255, 183, 3),
	sf::Color::Black,
	sf::Color(2, 48, 71)
	};
	bool paused = false;
};

#endif ViewH