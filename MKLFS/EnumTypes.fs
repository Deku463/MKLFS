module EnumTypes

type Layout =
    | RowMajor = 101
    | ColMajor = 102

type Trans =
    | NoTrans = 111
    | Trans = 112
    | ConjTrans = 113

type UpLo =
    | Upper = 121
    | Lower = 122

type Diag =
    | NonUnit = 131
    | Unit = 132

type Side =
    | Left = 141
    | Rigth = 142