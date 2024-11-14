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
	bool paused = false;
};

#endif ViewH