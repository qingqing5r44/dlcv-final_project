#include "main_window.h"
int main(int argc, char** argv) 
{
	QApplication a(argc, argv);
	Main* main_window = new Main;
	main_window->argc = argc;
	main_window->argv = argv;
	main_window->setAttribute(Qt::WA_DeleteOnClose);
	main_window->exec();
}