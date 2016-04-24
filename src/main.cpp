#include <QApplication>
#include <QDebug>

#include "maincontroller.h"
#include "networkstorage.h"
#include "undooperation.h"

int main(int argc, char** argv)
{
    NetworkStorage* st = new NetworkStorage("network");
    UndoOperation* undo = new UndoOperation("Undo");
    qDebug() << st->name();
    qDebug() << undo->name();
    QApplication app(argc, argv);

    int code = 0;
    {

        MainController controller;
        controller.show();
        code = app.exec();
    }
    delete undo;
    delete st;

    return code;
}
