#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include "Macierz.h"

#define ODCZYTY_GDANSK 512	
#define ODCZYTY_KANION 512
#define ODCZYTY_ULURU 512
struct odczyt {
	double x;
	double y;
};

void LU(Macierz* A, double* x, double* b);
void pivotLU(Macierz* U, Macierz* L, Macierz* P, int i);
double sklejanieFunkcji(odczyt* odczyty, double odleglosc, int iloscWezlow);
double interpolacjaLagrangea(odczyt* odczyty, double odleglosc, int iloscWezlow);
bool czytajCSV(const char* nazwaPliku, odczyt* odczyty, int ilosc);
bool czyZawiera(int* tab, int dlugosc, int wartosc);
int* losujWartosci(int max, int iloscWezlow);
double bladSredni(odczyt* oryginalne, odczyt* przyblizone, int ilosc);
void wykonajIZapisz(odczyt* odczyty, std::string nazwaPliku, int iloscWezlow, int iloscOdczytow, std::string metoda, std::string wybieranie);
bool czyWezel(odczyt* wezly, odczyt sprawdzany, int iloscWezlow);

int main()
{
	odczyt* odczyty = new odczyt[ODCZYTY_GDANSK];
	if(!czytajCSV("SpacerniakGdansk.csv", odczyty, ODCZYTY_GDANSK)) return 1;
	wykonajIZapisz(odczyty, "wyniki/GdanskLagrange5Wezlow.csv", 5, ODCZYTY_GDANSK, "lagrange", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/GdanskLagrange10Wezlow.csv", 10, ODCZYTY_GDANSK, "lagrange", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/GdanskLagrange15Wezlow.csv", 15, ODCZYTY_GDANSK, "lagrange", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/GdanskLagrange20Wezlow.csv", 20, ODCZYTY_GDANSK, "lagrange", "po kolei");

	wykonajIZapisz(odczyty, "wyniki/GdanskSklejanie5Wezlow.csv", 5, ODCZYTY_GDANSK, "sklejanie", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/GdanskSklejanie10Wezlow.csv", 10, ODCZYTY_GDANSK, "sklejanie", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/GdanskSklejanie15Wezlow.csv", 15, ODCZYTY_GDANSK, "sklejanie", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/GdanskSklejanie20Wezlow.csv", 20, ODCZYTY_GDANSK, "sklejanie", "po kolei");


	wykonajIZapisz(odczyty, "wyniki/GdanskLagrange5WezlowLos.csv", 5, ODCZYTY_GDANSK, "lagrange", "losowo");
	wykonajIZapisz(odczyty, "wyniki/GdanskLagrange10WezlowLos.csv", 10, ODCZYTY_GDANSK, "lagrange", "losowo");
	wykonajIZapisz(odczyty, "wyniki/GdanskLagrange15WezlowLos.csv", 25, ODCZYTY_GDANSK, "lagrange", "losowo");
	wykonajIZapisz(odczyty, "wyniki/GdanskLagrange20WezlowLos.csv", 20, ODCZYTY_GDANSK, "lagrange", "losowo");

	wykonajIZapisz(odczyty, "wyniki/GdanskSklejanie5WezlowLos.csv", 5, ODCZYTY_GDANSK, "sklejanie", "losowo");
	wykonajIZapisz(odczyty, "wyniki/GdanskSklejanie10WezlowLos.csv", 10, ODCZYTY_GDANSK, "sklejanie", "losowo");
	wykonajIZapisz(odczyty, "wyniki/GdanskSklejanie15WezlowLos.csv", 15, ODCZYTY_GDANSK, "sklejanie", "losowo");
	wykonajIZapisz(odczyty, "wyniki/GdanskSklejanie20WezlowLos.csv", 20, ODCZYTY_GDANSK, "sklejanie", "losowo");
	delete[] odczyty;

	odczyty = new odczyt[ODCZYTY_KANION];
	if(!czytajCSV("kolorado2.csv", odczyty, ODCZYTY_KANION)) return 1;
	wykonajIZapisz(odczyty, "wyniki/KanionLagrange5Wezlow.csv", 5, ODCZYTY_KANION, "lagrange", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/KanionLagrange10Wezlow.csv", 10, ODCZYTY_KANION, "lagrange", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/KanionLagrange15Wezlow.csv", 15, ODCZYTY_KANION, "lagrange", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/KanionLagrange20Wezlow.csv", 20, ODCZYTY_KANION, "lagrange", "po kolei");

	wykonajIZapisz(odczyty, "wyniki/KanionSklejanie5Wezlow.csv", 5, ODCZYTY_KANION, "sklejanie", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/KanionSklejanie10Wezlow.csv", 10, ODCZYTY_KANION, "sklejanie", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/KanionSklejanie15Wezlow.csv", 15, ODCZYTY_KANION, "sklejanie", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/KanionSklejanie20Wezlow.csv", 20, ODCZYTY_KANION, "sklejanie", "po kolei");


	wykonajIZapisz(odczyty, "wyniki/KanionLagrange5WezlowLos.csv", 5, ODCZYTY_KANION, "lagrange", "losowo");
	wykonajIZapisz(odczyty, "wyniki/KanionLagrange10WezlowLos.csv", 10, ODCZYTY_KANION, "lagrange", "losowo");
	wykonajIZapisz(odczyty, "wyniki/KanionLagrange15WezlowLos.csv", 15, ODCZYTY_KANION, "lagrange", "losowo");
	wykonajIZapisz(odczyty, "wyniki/KanionLagrange20WezlowLos.csv", 20, ODCZYTY_KANION, "lagrange", "losowo");

	wykonajIZapisz(odczyty, "wyniki/KanionSklejanie5WezlowLos.csv", 5, ODCZYTY_KANION, "sklejanie", "losowo");
	wykonajIZapisz(odczyty, "wyniki/KanionSklejanie10WezlowLos.csv", 10, ODCZYTY_KANION, "sklejanie", "losowo");
	wykonajIZapisz(odczyty, "wyniki/KanionSklejanie15WezlowLos.csv", 15, ODCZYTY_KANION, "sklejanie", "losowo");
	wykonajIZapisz(odczyty, "wyniki/KanionSklejanie20WezlowLos.csv", 20, ODCZYTY_KANION, "sklejanie", "losowo");
	delete[] odczyty;

	odczyty = new odczyt[ODCZYTY_ULURU];
	if(!czytajCSV("uluru.csv", odczyty, ODCZYTY_ULURU)) return 1;
	wykonajIZapisz(odczyty, "wyniki/UluruLagrange5Wezlow.csv", 5, ODCZYTY_ULURU, "lagrange", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/UluruLagrange10Wezlow.csv", 10, ODCZYTY_ULURU, "lagrange", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/UluruLagrange15Wezlow.csv", 15, ODCZYTY_ULURU, "lagrange", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/UluruLagrange20Wezlow.csv", 20, ODCZYTY_ULURU, "lagrange", "po kolei");

	wykonajIZapisz(odczyty, "wyniki/UluruSklejanie5Wezlow.csv", 5, ODCZYTY_ULURU, "sklejanie", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/UluruSklejanie10Wezlow.csv", 10, ODCZYTY_ULURU, "sklejanie", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/UluruSklejanie15Wezlow.csv", 15, ODCZYTY_ULURU, "sklejanie", "po kolei");
	wykonajIZapisz(odczyty, "wyniki/UluruSklejanie20Wezlow.csv", 20, ODCZYTY_ULURU, "sklejanie", "po kolei");


	wykonajIZapisz(odczyty, "wyniki/UluruLagrange5WezlowLos.csv", 5, ODCZYTY_ULURU, "lagrange", "losowo");
	wykonajIZapisz(odczyty, "wyniki/UluruLagrange10WezlowLos.csv", 10, ODCZYTY_ULURU, "lagrange", "losowo");
	wykonajIZapisz(odczyty, "wyniki/UluruLagrange15WezlowLos.csv", 15, ODCZYTY_ULURU, "lagrange", "losowo");
	wykonajIZapisz(odczyty, "wyniki/UluruLagrange20WezlowLos.csv", 20, ODCZYTY_ULURU, "lagrange", "losowo");

	wykonajIZapisz(odczyty, "wyniki/UluruSklejanie5WezlowLos.csv", 5, ODCZYTY_ULURU, "sklejanie", "losowo");
	wykonajIZapisz(odczyty, "wyniki/UluruSklejanie10WezlowLos.csv", 10, ODCZYTY_ULURU, "sklejanie", "losowo");
	wykonajIZapisz(odczyty, "wyniki/UluruSklejanie15WezlowLos.csv", 15, ODCZYTY_ULURU, "sklejanie", "losowo");
	wykonajIZapisz(odczyty, "wyniki/UluruSklejanie20WezlowLos.csv", 20, ODCZYTY_ULURU, "sklejanie", "losowo");

	delete[] odczyty;
}


bool czyZawiera(int* tab, int dlugosc, int wartosc) {
	for (int i = 0; i < dlugosc; i++) {
		if (tab[i] == wartosc) return true;
	}
	return false;
}

int* losujWartosci(int max, int iloscWezlow) {
	int* tab = new int[iloscWezlow];

	tab[0] = 0;
	for (int i = 1; i < iloscWezlow - 1; i++) {
		int r = rand() % (max - 2) + 1;
		if (!czyZawiera(tab, iloscWezlow, r))tab[i] = r;
		else i--;
	}
	tab[iloscWezlow - 1] = max - 1;
	for (int i = 0; i < iloscWezlow - 1; i++) {

		for (int j = 0; j < iloscWezlow - i - 1; j++) {
			if (tab[j] > tab[j + 1])
				std::swap(tab[j], tab[j + 1]);
		}
	}

	return tab;
}

double bladSredni(odczyt* oryginalne, odczyt* przyblizone, int ilosc) {

	double blad = 0;
	for (int i = 0; i < ilosc; i++) {
		blad += abs(oryginalne[i].y - przyblizone[i].y);
	}
	return blad / ilosc;
}

void wykonajIZapisz(odczyt* odczyty, std::string nazwaPliku, int iloscWezlow, int iloscOdczytow, std::string metoda, std::string wybieranie) {
	auto start = std::chrono::high_resolution_clock::now();
	odczyt* wybraneWartosci = new odczyt[iloscWezlow];

	if (wybieranie == "po kolei") {
		int krok = (iloscOdczytow - 2) / (iloscWezlow - 1);
		int nr = 0;
		for (int i = 0; i < iloscWezlow - 1; i++) {

			wybraneWartosci[i] = odczyty[nr];
			nr += krok;
		}
		wybraneWartosci[iloscWezlow - 1] = odczyty[iloscOdczytow - 1];
	}
	else if (wybieranie == "losowo") {
		int* losowe = losujWartosci(iloscOdczytow, iloscWezlow);
		for (int i = 0; i < iloscWezlow; i++) {
			wybraneWartosci[i] = odczyty[losowe[i]];
		}
		delete[] losowe;
	}
	std::ofstream file;
	file.open(nazwaPliku);

	file << "Dystans (m),Wysokosc (m),Wezel\n";
	double blad = 0;
	if (metoda == "lagrange") {
		for (int i = 0; i < iloscOdczytow; i++) {
			double ilg = interpolacjaLagrangea(wybraneWartosci, odczyty[i].x, iloscWezlow);
			bool isNode = czyWezel(wybraneWartosci, odczyty[i], iloscWezlow);
			file << odczyty[i].x << "," << ilg << ","<< isNode <<"\n";
			blad += abs(ilg - odczyty[i].y);
		}

	}
	else if (metoda == "sklejanie") {
		for (int i = 0; i < iloscOdczytow; i++) {
			double sf = sklejanieFunkcji(wybraneWartosci, odczyty[i].x, iloscWezlow);
			bool isNode = czyWezel(wybraneWartosci, odczyty[i], iloscWezlow);
			file << odczyty[i].x << "," << sf << "," << isNode << "\n";
			blad += abs(sf - odczyty[i].y);
		}
	}
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	blad /= iloscOdczytow;
	std::cout << nazwaPliku << " - blad sredni: " << blad << ", czas: " << elapsed.count() << " s\n";
	delete[] wybraneWartosci;
	file.close();
}

double interpolacjaLagrangea(odczyt* odczyty, double odleglosc, int iloscWezlow) {

	double y = 0;
	for (int i = 0; i < iloscWezlow; i++) {
		double iloczyn = 1;
		for (int j = 0; j < iloscWezlow; j++) {
			if (i != j) iloczyn *= (odleglosc - odczyty[j].x) / (odczyty[i].x - odczyty[j].x);
		}
		y += iloczyn * odczyty[i].y;
	}

	return y;
}

bool czytajCSV(const char* nazwaPliku, odczyt* odczyty, int ilosc) {
	std::ifstream file;
	file.open(nazwaPliku, file.in);

	if (!file.is_open()) {
		printf_s("Niepowodzenie przy otwieraniu pliku!");
		return false;
	}

	std::string dane;
	getline(file, dane);
	for (int i = 0; i < ilosc; i++) {
		getline(file, dane, ',');
		odczyty[i].x = stod(dane);
		getline(file, dane);
		odczyty[i].y = stod(dane);
	}

	file.close();
	return true;
}

void pivotLU(Macierz* U, Macierz* L, Macierz* P, int i) {
	double pivot = abs(U->tab[i][i]);
	int index = i;

	for (int j = i + 1; j < U->rozmiar; j++) {
		if (abs(U->tab[j][i]) > pivot) {
			pivot = abs(U->tab[j][i]);
			index = j;
		}
	}

	if (index != i) {
		double pom;

		for (int j = 0; j < U->rozmiar; j++)	{
			if (j >= i)	{
				pom = U->tab[i][j];
				U->tab[i][j] = U->tab[index][j];
				U->tab[index][j] = pom;
			}
			else {
				pom = L->tab[i][j];
				L->tab[i][j] = L->tab[index][j];
				L->tab[index][j] = pom;
			}

			pom = P->tab[i][j];
			P->tab[i][j] = P->tab[index][j];
			P->tab[index][j] = pom;
		}
	}
}

void LU(Macierz* A, double* x, double* b) {
	Macierz* U = A->kopia();
	Macierz* L = new Macierz(A->rozmiar, 1, 0, 0);
	Macierz* P = new Macierz(A->rozmiar, 1, 0, 0);

	for (int i = 0; i < A->rozmiar - 1; i++) {

		pivotLU(U, L, P, i);

		for (int j = i + 1; j < A->rozmiar; j++) {
			L->tab[j][i] = U->tab[j][i] / U->tab[i][i];

			for (int k = i; k < A->rozmiar; k++) {
				U->tab[j][k] -= L->tab[j][i] * U->tab[i][k];
			}
		}
	}

	double* y = new double[A->rozmiar];
	for (int i = 0; i < A->rozmiar; i++) {
		double suma = 0;
		for (int j = 0; j < i; j++) {
			suma += L->tab[i][j] * y[j];
		}
		y[i] = (b[i] - suma) / L->tab[i][i];
	}

	for (int i = A->rozmiar - 1; i >= 0; i--) {
		double suma = 0;
		for (int j = i + 1; j < A->rozmiar; j++) {
			suma += U->tab[i][j] * x[j];
		}
		x[i] = (y[i] - suma) / U->tab[i][i];
	}
	delete U;
	delete L;
	delete P;
	delete[] y;
}

double sklejanieFunkcji(odczyt* odczyty, double dystans, int iloscWezlow) {

	Macierz* M = new Macierz(4 * (iloscWezlow - 1));
	double* x = new double[4 * (iloscWezlow - 1)];
	double* b = new double[4 * (iloscWezlow - 1)];

	double h = odczyty[1].x - odczyty[0].x;
	M->tab[0][0] = 1;
	b[0] = odczyty[0].y;
	M->tab[1][0] = 1;
	M->tab[1][1] = h;
	M->tab[1][2] = pow(h, 2);
	M->tab[1][3] = pow(h, 3);
	b[1] = odczyty[1].y;
	M->tab[2][2] = 1;
	b[2] = 0;
	h = odczyty[iloscWezlow - 1].x - odczyty[iloscWezlow - 2].x;
	M->tab[3][4 * (iloscWezlow - 2) + 2] = 2;
	M->tab[3][4 * (iloscWezlow - 2) + 3] = 6 * h;
	b[3] = 0;

	for (int i = 1; i < iloscWezlow - 1; i++) {
		M->tab[4 * i][4 * i] = 1;
		b[4 * i] = odczyty[i].y;

		h = odczyty[i + 1].x - odczyty[i].x;
		M->tab[4 * i + 1][4 * i] = 1;
		M->tab[4 * i + 1][4 * i + 1] = h;
		M->tab[4 * i + 1][4 * i + 2] = pow(h, 2);
		M->tab[4 * i + 1][4 * i + 3] = pow(h, 3);
		b[4 * i + 1] = odczyty[i + 1].y;

		M->tab[4 * i + 2][4 * (i - 1) + 1] = 1;
		M->tab[4 * i + 2][4 * (i - 1) + 2] = 2 * h;
		M->tab[4 * i + 2][4 * (i - 1) + 3] = 3 * pow(h, 2);
		M->tab[4 * i + 2][4 * i + 1] = -1;
		b[4 * i + 2] = 0;

		M->tab[4 * i + 3][4 * (i - 1) + 2] = 2;
		M->tab[4 * i + 3][4 * (i - 1) + 3] = 6 * h;
		M->tab[4 * i + 3][4 * i + 2] = -2;
		b[4 * i + 3] = 0;
	}
	LU(M, x, b);

	double y = 0;
	for (int i = 0; i < iloscWezlow - 1; i++) {
		if ((dystans >= odczyty[i].x) && (dystans <= odczyty[i + 1].x)) {
			for (int j = 0; j < 4; j++) {
				h = dystans - odczyty[i].x;
				y += x[4 * i + j] * pow(h, j);
			}
			break;
		}
	}
	delete[] x;
	delete[] b;
	delete M;

	return y;
}

bool czyWezel(odczyt* wezly, odczyt sprawdzany, int iloscWezlow) {
	for (int i = 0; i < iloscWezlow; i++) {
		if (sprawdzany.x == wezly[i].x) return true;
	}
	return false;
}