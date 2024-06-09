#include "Batchview.h"

BatchView::BatchView(sf::RenderWindow& window, Landscape* pLand, Community* pCommunity) :
	pLandscape{ pLand },
	pComm{ pCommunity } {

	// Set cell size such that both dimensions fit on screen
	dimX = pLandscape->getLandParams().dimX;
	dimY = pLandscape->getLandParams().dimY;
	cellSize = min(1920u / dimX, 1080u / dimY);

	unsigned int winWidth = dimX * cellSize;
	unsigned int winHeight = dimY * cellSize + dimY * relSizeLegend; // space at bottom for legend

	window.create(sf::VideoMode{winWidth, winHeight}, "RangeShifter Batch");
	window.setFramerateLimit(144);
}


// ---------------------------------------
// Collect and process user window input
// e.g. clicking, scrolling, typing etc.
// ---------------------------------------
void BatchView::collectUserInput(sf::RenderWindow& window) {
	if (window.isOpen()) {
		for (auto event = sf::Event{}; window.pollEvent(event);)
		{
			// Close the window
			if (event.type == sf::Event::Closed)
			{
				window.close();
			}
			// else, add more events...
		}
	}
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
			if (pCell != 0) {
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
	for (int iPch : patchIndices) {
		if (iPch == 0) continue; // ignore individuals in matrix
		pPatch = pLandscape->findPatch(iPch);
		pPop = (Population*)pPatch->getPopn((intptr)pSpecies);
		if (pPop != 0) {
			popSize = pPop->getNInds();
			for (int i = 0; i < popSize; i++) {
				pRandCell = pPatch->getRandomCell();
				randLocn = pRandCell->getLocn();

				cInd = indShape;
				cInd.setFillColor(sf::Color::Blue);

				if (pPatch->getPatchNum() == 30) {
					cInd.setFillColor(indColour);

					patchLimits plim = pPatch->getLimits();
					sf::RectangleShape thatPatch{sf::Vector2f(plim.xMax - plim.xMin, plim.yMax - plim.yMin)};
					thatPatch.setOutlineColor(indColour);
					thatPatch.setOutlineThickness(3.0);
					thatPatch.setFillColor(sf::Color::Transparent);
					thatPatch.setPosition(plim.xMin, plim.xMax);
					window.draw(thatPatch);
				}
				
				// Randomise position inside the cell
				indPosition = { (float)((randLocn.x * pRandom->Random() * cellSize)),
					float((randLocn.y + pRandom->Random()) * cellSize) };
				cInd.setPosition(indPosition);
				window.draw(cInd);
			}
		}
	}

	// Display year and generation
	sf::Vector2f timeLegendPosition{0.0, (1 + relSizeLegend) * dimY}; // just below plot
	sf::RectangleShape timeLegendBg(sf::Vector2f{(float)dimX, relSizeLegend * dimY});
	timeLegendBg.setPosition(timeLegendPosition);
	sf::Text timeLegendTxt;
	sf::Font font;
	font.loadFromFile("../../../graphic_resources/rockwell.otf");
	timeLegendTxt.setFont(font);
	timeLegendTxt.setPosition(timeLegendPosition);
	timeLegendTxt.setString("Year:\t" + to_string(yr) + "\n Gen:\t" + to_string(gen));
	timeLegendTxt.setFillColor(sf::Color::White);
	window.draw(timeLegendBg);
	window.draw(timeLegendTxt);

	window.display();
}