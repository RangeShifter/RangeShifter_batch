#include "Batchview.h"
#include <filesystem>

BatchView::BatchView(sf::RenderWindow& window, Landscape* pLand, Community* pCommunity,
	const int& nbYears, const int& nbSeasons) :
	pLandscape{ pLand },
	pComm{ pCommunity },
	maxYear { nbYears },
	maxSeason { nbSeasons }
{
	string pathToFont = "../../../gfx/consola.ttf";
	if (!std::filesystem::exists(pathToFont)) {
		throw logic_error("Font doesn't exist.\n");

	}
	if (!font.loadFromFile(pathToFont)) {
		throw logic_error("Font couldn't be loaded.\n");
	}

	// Set cell size such that both dimensions fit on screen
	dimX = pLandscape->getLandParams().dimX;
	dimY = pLandscape->getLandParams().dimY;
	cellSize = min(1920u / dimX, 1080u / dimY);

	unsigned int winWidth = dimX * cellSize;
	unsigned int winHeight = dimY * cellSize 
		+ dimY * relSizeLegend; // space at bottom for legend

	window.create(sf::VideoMode{winWidth, winHeight}, "RangeShifter Batch");
	window.setFramerateLimit(144);

	// Create legend info text box
	float legendOriginY = (1.0 - relSizeLegend) * dimY;
	float legendOriginX = 0.01f * dimX;
	sf::Vector2f timeLegendPosition{ 0.0, legendOriginY }; // just below plot
	timeLegendBg = sf::RectangleShape(sf::Vector2f{ dimX * 1.0f, legendOriginY });
	timeLegendBg.setPosition(timeLegendPosition);

	timeLegendTxt.setFont(font);
	timeLegendTxt.setPosition(sf::Vector2f(legendOriginX, legendOriginY));
	timeLegendTxt.setCharacterSize(relSizeLegend * dimY / 2.0f);
	timeLegendTxt.setColor(sf::Color::Blue);

	// Create paused text box
	float txtPausedX = dimX / 2.0f;
	float txtPausedY = dimY / 2.0f;
	txtPaused.setFont(font);
	txtPaused.setOrigin(0.5, 0.5);
	txtPaused.setPosition(sf::Vector2f(txtPausedX, txtPausedY));
	txtPaused.setString("[Paused]");
	txtPaused.setColor(sf::Color::Blue);
	txtPausedBg = sf::RectangleShape(sf::Vector2f{
		txtPaused.getGlobalBounds().width * 1.1f,
		txtPaused.getGlobalBounds().height * 1.4f
		});
	txtPausedBg.setPosition(sf::Vector2f(txtPausedX, txtPausedY));
	txtPausedBg.setFillColor(sf::Color::White);
}


// ---------------------------------------
// Collect and process user window input
// e.g. clicking, scrolling, typing etc.
// ---------------------------------------
void BatchView::collectUserInput(sf::RenderWindow& window) {

	for (auto event = sf::Event{}; window.pollEvent(event);) {
		switch (event.type) {
			// Close the window
		case sf::Event::Closed:
			window.close();
			break;
			// Pause the simulation
		case sf::Event::KeyPressed:
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Space)) {
				paused = !paused; // switch on/off
			}
			break;
		default:
			break;
		}
	}
} 

void BatchView::runPauseLoop(sf::RenderWindow& window) {
	while (paused && window.isOpen()) {
		collectUserInput(window);
	};
}

void BatchView::drawLandscape(sf::RenderWindow& window) {
	const int maxY = pLandscape->getLandParams().maxY * cellSize;

	for (int x = 0; x < dimX; x++) {

		// Displaying can take a while, 
		// check input between columns
		collectUserInput(window);
		if (!window.isOpen()) return; // dinnae bother finishing

		for (int y = 0; y < dimY; y++) {

			Cell* pCell = pLandscape->findCell(x, y);
			if (pCell != nullptr) {
				sf::RectangleShape c(sf::Vector2f(cellSize, cellSize));
				int h = pCell->getHabIndex(0);
				sf::Color col = habitatPalette[h];
				c.setPosition(cellSize * x, maxY - cellSize * y);
				c.setFillColor(col);
				window.draw(c);
			}
		}
	}

}

void BatchView::drawCommunity(sf::RenderWindow& window, Species* pSpecies, const int& yr, const int& gen) {

	// Erase previous community
	window.clear();
	drawLandscape(window);

	// Draw current community
	Patch* pPatch;
	Population* pPop;
	int popSize;
	Cell* pRandCell;
	locn randLocn;
	sf::Vector2f indPosition;

	sf::CircleShape cInd;
	
	const vector<int> patchIndices = pLandscape->readPatchNums();
	patchLimits plim;

	for (int iPch : patchIndices) {

		if (iPch == 0) continue; // ignore individuals in matrix
		pPatch = pLandscape->findPatch(iPch);
		pPop = (Population*)pPatch->getPopn((intptr)pSpecies);

		if (pPop != nullptr) {
			popSize = pPop->getNInds();

			for (int i = 0; i < popSize; i++) {
				pRandCell = pPatch->getRandomCell();
				randLocn = pRandCell->getLocn();

				cInd = indShape;
				cInd.setFillColor(sf::Color::Blue);

				if (pPatch->getPatchNum() == 30) {

					cInd.setFillColor(indColour);

					plim = pPatch->getLimits();
					sf::RectangleShape thatPatch{sf::Vector2f(plim.xMax - plim.xMin, plim.yMax - plim.yMin)};
					thatPatch.setOutlineColor(indColour);
					thatPatch.setOutlineThickness(3.0);
					thatPatch.setFillColor(sf::Color::Transparent);
					thatPatch.setPosition(plim.xMin, plim.xMax);
					//window.draw(thatPatch);
				}
				
				// Randomise position inside the cell
				indPosition = { 
					static_cast<float>(randLocn.x * pRandom->Random() * cellSize),
					static_cast<float>(randLocn.y + pRandom->Random() * cellSize)
				};
				cInd.setPosition(indPosition);
				window.draw(cInd);
			}
		}
	}

	// Display year and season
	timeLegendTxt.setString(
		"Year:\t" + to_string(yr + 1) + "/" + to_string(maxYear) +
		+"\nSeason:\t" + to_string(gen + 1) + "/" + to_string(maxSeason)
	);

	window.draw(timeLegendBg);
	window.draw(timeLegendTxt);

	if (paused || yr + 1 == maxYear) {
		window.draw(txtPausedBg);
		window.draw(txtPaused);
	}

	window.display();
}