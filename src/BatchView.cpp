#include "Batchview.h"

BatchView::BatchView(sf::RenderWindow& window, Landscape* pLand, Community* pCommunity,
	const int& nbYears, const int& nbSeasons) :
	pLandscape{ pLand },
	pComm{ pCommunity },
	maxYear { nbYears },
	maxSeason { nbSeasons }
{
	const string pathToFont = "gfx/consola.ttf";
	if (!std::filesystem::exists(pathToFont)) {
		throw logic_error("Font file is missing. Please copy " 
			+ pathToFont + " into the project directory.\n");

	} else if (!font.loadFromFile(pathToFont)) {
		throw logic_error("Font file exists but could not be loaded.\n");
	}

	// Set cell size such that both dimensions fit on screen
	const int dimX = pLandscape->getLandParams().dimX;
	const int dimY = pLandscape->getLandParams().dimY;
	cellSize = min( 
		(1.0f * dfltWinWidth) / dimX, 
		dfltWinHeight * (1 - relSizeLegend) / dimY
	);

	winWidth = dimX * cellSize;
	winHeight = dimY * cellSize 
		/ (1 - relSizeLegend); // space at bottom for legend

	window.create(sf::VideoMode{winWidth, winHeight}, "RangeShifter Batch");
	window.setFramerateLimit(144);

	// Create legend info text box
	float legendOriginY = (1.0 - relSizeLegend) * winHeight;
	float legendOriginX = 0.01f * dimX;
	sf::Vector2f timeLegendPosition{ 0.0, legendOriginY }; // just below plot
	timeLegendBg = sf::RectangleShape(sf::Vector2f{ dimX * 1.0f, legendOriginY });
	timeLegendBg.setPosition(timeLegendPosition);

	timeLegendTxt.setFont(font);
	timeLegendTxt.setPosition(sf::Vector2f(legendOriginX, legendOriginY));
	timeLegendTxt.setCharacterSize(winHeight * relSizeLegend / 3.0f);
	timeLegendTxt.setColor(sf::Color::Blue);

	// Create paused text box
	txtPaused.setFont(font);
	txtPaused.setString("[Paused]");
	float txtPausedX = winWidth / 2.0f - txtPaused.getGlobalBounds().width / 2.0f;
	float txtPausedY = 0.0; // legendOriginY;
	txtPaused.setPosition(sf::Vector2f(txtPausedX, txtPausedY));
	txtPaused.setColor(sf::Color::Blue);
	float bgWidth = txtPaused.getGlobalBounds().width * 1.2f;
	float bgHeight = txtPaused.getGlobalBounds().height * 1.4f;
	txtPausedBg = sf::RectangleShape(sf::Vector2f{bgWidth, bgHeight});
	txtPausedBg.setPosition(sf::Vector2f(
		winWidth / 2.0f - bgWidth / 2.0f,
		txtPausedY
	));
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
	const int dimX = pLandscape->getLandParams().dimX;
	const int dimY = pLandscape->getLandParams().dimY;

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
	Population* pPop = nullptr;
	sf::CircleShape cInd;
	
	const vector<int> patchIndices = pLandscape->readPatchNums();
	for (int iPch : patchIndices) {

		if (iPch == 0) continue; // ignore individuals in matrix
		pPatch = pLandscape->findPatch(iPch);
		pPop = (Population*)pPatch->getPopn((intptr)pSpecies);

		if (pPop != nullptr) {

			//vector<locn> juvCoords = pPop->getIndsCoords(true);
			vector<locn> indCoords = pPop->getIndsCoords(false);

			/*
			for (auto coord : juvCoords) {
				cInd = indShape;
				cInd.setFillColor();
				cInd.setPosition(coord.x * cellSize * 1.0f, coord.y * cellSize * 1.0f);
				window.draw(cInd);
			}
			*/
			for (auto coord : indCoords) {
				cInd = indShape;
				cInd.setFillColor(indColour);
				cInd.setPosition(coord.x * cellSize * 1.0f, coord.y * cellSize * 1.0f);
				window.draw(cInd);
			} 
		}
	}

	// Display year and season
	timeLegendTxt.setString(
		"Year: " + to_string(yr + 1) + "/" + to_string(maxYear) +
		+"\nSeason: " + to_string(gen + 1) + "/" + to_string(maxSeason)
	);

	window.draw(timeLegendBg);
	window.draw(timeLegendTxt);

	if (paused || yr + 1 == maxYear) {
		window.draw(txtPausedBg);
		window.draw(txtPaused);
	}

	window.display();
}