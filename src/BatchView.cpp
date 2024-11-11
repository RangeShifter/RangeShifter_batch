#include "Batchview.h"
#include <filesystem>

BatchView::BatchView(sf::RenderWindow& window, Landscape* pLand, Community* pCommunity) :
	pLandscape{ pLand },
	pComm{ pCommunity } 
{

	string pathToFont = "../../../gfx/arial.ttf";

	cout << std::filesystem::current_path() << endl;

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

}


// ---------------------------------------
// Collect and process user window input
// e.g. clicking, scrolling, typing etc.
// ---------------------------------------
void BatchView::collectUserInput(sf::RenderWindow& window) {
	mustPause = false;
	do {
		if (window.isOpen()) {
			for (auto event = sf::Event{}; window.pollEvent(event);) {
				switch (event.type) {
				// Close the window
				case sf::Event::Closed:
					window.close();
					break;
				// Pause the simulation
				case sf::Event::KeyPressed:
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::Space)) {
						mustPause = !mustPause; // pause if not paused and conversely
						if (mustPause) {
							//sf::Text pauseMsg;
						}
					}
					break;
				default:
					break;
				}
			}
		}
		else break; // ignore pause if window has closed
	} while (mustPause);
	
}

void BatchView::drawLandscape(sf::RenderWindow& window) {
	const int maxY = pLandscape->getLandParams().maxY * cellSize;

	for (int x = 0; x < dimX; x++) {

		// Displaying can take long, 
		// check keys/buttons between columns
		collectUserInput(window);
		if (!window.isOpen()) return; // dinnae bother

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
					window.draw(thatPatch);
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

	// Display year and generation
	float legendHeight = (1.0 - relSizeLegend) * dimY;
	sf::Vector2f timeLegendPosition{0.0, legendHeight}; // just below plot
	sf::RectangleShape timeLegendBg(sf::Vector2f{dimX * 1.0f, legendHeight});
	timeLegendBg.setPosition(timeLegendPosition);
	window.draw(timeLegendBg);

	sf::Text timeLegendTxt;
	timeLegendTxt.setFont(font);
	timeLegendTxt.setPosition(timeLegendPosition);
	timeLegendTxt.setString(
		"Year:\t" + to_string(yr) 
		+ "\n Gen:\t" + to_string(gen)
	);
	timeLegendTxt.setColor(sf::Color::Blue);
	window.draw(timeLegendTxt);

	sf::Text hello;
	hello.setFont(font);
	hello.setPosition(sf::Vector2f{234.0, 353.0 });
	hello.setCharacterSize(100);
	hello.setString("Hello");
	hello.setColor(sf::Color::Blue);
	window.draw(hello);

	window.display();
}