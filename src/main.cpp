#include <QApplication>
#include <QDebug>

#include "maincontroller.h"
#include "networkstorage.h"

int main(int argc, char** argv)
{
    NetworkStorage* st = new NetworkStorage("network");
    qDebug() << st->name();
    QApplication app(argc, argv);

    int code = 0;
    {

        MainController controller;
        controller.show();
        code = app.exec();
    }
    delete st;
    return code;
}
