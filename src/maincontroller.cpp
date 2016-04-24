#include "maincontroller.h"

MainController::MainController()
    :m_window(NULL)
{

}

MainController::~MainController()
{

}

void MainController::show()
{
    m_window = new MainWindow;
    m_window->show();
}
