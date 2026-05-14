**Что Такое ProtT5**  
ProtT5 - это protein language model из проекта **ProtTrans** от Rostlab/TUM. Он основан на архитектуре [[T5]], обучен на белковых последовательностях, а не на человеческом тексте.

ProtT5-XL-UniRef50:
- около 3B параметров в полной версии;
- обычно используют encoder-only вариант;
- выдает [[Embedding]] L x 1024, где L - длина белка;
- потом часто делают mean pooling и получают один 1024-d вектор на белок.

Модель есть на Hugging Face: [Rostlab/prot_t5_xl_half_uniref50-enc](https://huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc).