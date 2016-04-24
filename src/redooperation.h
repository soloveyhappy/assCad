#ifndef REDOOPERATION_H
#define REDOOPERATION_H

#include "operation.h"

class RedoOperation : public Operation
{
public:
    using Operation::Operation;

    // Operation interface
public:
    bool cancel() override;
    bool execute() override;
};

#endif // REDOOPERATION_H
