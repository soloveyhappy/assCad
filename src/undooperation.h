#ifndef UNDOOPERATION_H
#define UNDOOPERATION_H

#include "operation.h"


class UndoOperation : public Operation
{
public:
    using Operation::Operation;
public:
    bool cancel() override;
    bool execute() override;
};

#endif // UNDOOPERATION_H
