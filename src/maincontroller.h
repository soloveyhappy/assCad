#ifndef MAINCONTROLLER_H
#define MAINCONTROLLER_H

#include <QObject>
#include "mainwindow.h"

class MainController : public QObject
{
    Q_OBJECT
public:
    MainController();
    ~MainController();
public:
    void show();

private:
    MainWindow* m_window;

signals:

public slots:
};

#endif // MAINCONTROLLER_H
