#include "Batchview.h"

BatchView::BatchView(sf::RenderWindow& window, Landscape* pLand, Community* pCommunity) :
	pLandscape{ pLand },
	pComm{ pCommunity } {

	// Set cell size such that both dimensions fit on screen
	dimX = pLandscape->getLandParams().dimX;
	dimY = pLandscape->getLandParams().dimY;
	cellSize = min(1920u / dimX, 1080u / dimY);

	window.create(sf::VideoMode{ dimX * cellSize, dimY * cellSize }, "RangeShifter Batch");
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

void BatchView::drawCommunity(sf::RenderWindow& window) {

	// Erase previous community
	window.clear();
	drawLandscape(window);

	// Draw current community
	// ...
	// for (each patch)
	// getPopn()
	// for each ind in Popn
	// get random cell in Patch
	// get random x and y in cell
	// draw individual shape there

	window.display();
}