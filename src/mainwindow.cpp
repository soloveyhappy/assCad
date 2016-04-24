#include "mainwindow.h"
#include <QWidget>
#include <QDockWidget>
#include <QStatusBar>
#include <QLayout>

MainWindow::MainWindow()
{

    QWidget* mainWidget = new QWidget ;
    QStatusBar* status = new QStatusBar;
    QDockWidget dataDoc;

    QVBoxLayout vertLayout;


    setStatusBar(status);


    setCentralWidget(mainWidget);
}

MainWindow::~MainWindow()
{

}
