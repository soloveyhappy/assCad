#include "mainwindow.h"
#include <QWidget>
#include <QDockWidget>
#include <QStatusBar>
#include <QLayout>
#include <QMouseEvent>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QApplication>
#include <QDockWidget>
#include <QTreeView>
#include <QTextEdit>
#include <QTreeView>
#include <QDebug>
#include <QSplitter>

#include "treemodel.h"

MainWindow::MainWindow()
{

    //QWidget* mainWidget = new QWidget ;
    text = new QTextEdit;
    QString s = "Getting Started  \t How to familiarize yourself with Qt Designer\n"
                "   Launching Designer   \t   Running the Qt Designer application\n"
                "The User Interface  \t   How to interact with Qt Designer\n"
                "   The User Interface   \t   How to interact with Qt Designer\n";

    text->setPlainText(s);


    QPushButton* btn = new QPushButton("sdsd");
    connect(btn, SIGNAL(clicked(bool)), this, SLOT(clicked(bool)));
    QVBoxLayout* layout = new QVBoxLayout;
    layout->addWidget(text);
    layout->addWidget(btn);
    QWidget* centralWidget = new QWidget;

    centralWidget->setLayout(layout);
    setCentralWidget(centralWidget);
    dataDoc = new QDockWidget( tr("Data"));
    addDockWidget(Qt::DockWidgetArea::LeftDockWidgetArea, dataDoc);
    setFixedSize(800, 600);


  //  QStatusBar* status = new QStatusBar;

//    QTreeView* dataTreeView = new QTreeView;



//    QVBoxLayout* vertLayout = new QVBoxLayout;

//    text = new QLabel("text");
//    vertLayout->addWidget(text);

//    QLabel* statusLabel = new QLabel("empty");
//    status->addPermanentWidget(statusLabel);
//    setStatusBar(status);

//
//
//    vertLayout->addWidget(btn);

    //mainWidget->setLayout(vertLayout);


}

MainWindow::~MainWindow()
{

}

void MainWindow::clicked(bool)
{
    QApplication::sendEvent(this, new QMouseEvent(QEvent::MouseButtonPress,
                                                  QPointF(10,10),
                                                  Qt::MouseButton::LeftButton,
                                                  Qt::MouseButton::LeftButton,
                                                  Qt::KeyboardModifier::NoModifier ));
    QString s = text->toPlainText();
    TreeModel* model = new TreeModel(s);
    QTreeView* modelView = new QTreeView;
    modelView->setModel(model);
    for (int column = 0; column < model->columnCount(); ++column)
            modelView->resizeColumnToContents(column);
    modelView->resizeColumnToContents(0);
    modelView->resizeColumnToContents(1);


    dataDoc->setWidget(modelView);
    dataDoc->resize(modelView->size());

}

void MainWindow::mousePressEvent(QMouseEvent* mouseEvent)
{
    QString s = QString("x:%1 - y:%2").arg(mouseEvent->localPos().x()).arg(mouseEvent->localPos().y());
    qDebug() << s;
    mouseEvent->setAccepted(false);
}

