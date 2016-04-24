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
    QString s = "Getting Started                         How to familiarize yourself with Qt Designer\
                Launching Designer                  Running the Qt Designer application \
                The User Interface                  How to interact with Qt Designer";
    text->setText(s);


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
    dataDoc->setWidget(modelView);

}

void MainWindow::mousePressEvent(QMouseEvent* mouseEvent)
{
    QString s = QString("x:%1 - y:%2").arg(mouseEvent->localPos().x()).arg(mouseEvent->localPos().y());
    qDebug() << s;
    mouseEvent->setAccepted(true);
}

