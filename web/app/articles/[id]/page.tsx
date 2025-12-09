'use client';

import { useEffect, useState } from 'react';
import { useParams, useRouter } from 'next/navigation';
import { api, ArticleViewModel } from '@/lib/api';
import StatusBadge from '@/components/StatusBadge';
import TabContainer from '@/components/TabContainer';
import ArticleViewer from '@/components/ArticleViewer';
import HumanFeedbackPanel from '@/components/HumanFeedbackPanel';
import LoadingSpinner from '@/components/LoadingSpinner';
import { useArticleStream } from '@/hooks/useArticleStream';
import { format } from 'date-fns';
import { ja } from 'date-fns/locale/ja';

export default function ArticleDetailPage() {
  const params = useParams();
  const router = useRouter();
  const articleId = params.id as string;

  const [article, setArticle] = useState<ArticleViewModel | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [executing, setExecuting] = useState<string | null>(null);
  const [successMessage, setSuccessMessage] = useState<string | null>(null);
  
  // ストリーミングフック
  const { events, connected } = useArticleStream(articleId);

  useEffect(() => {
    loadArticle();
  }, [articleId]);
  
  // ストリーミングイベントで記事を更新
  useEffect(() => {
    if (events.length > 0) {
      const lastEvent = events[events.length - 1];
      if (lastEvent.type === 'phase_change') {
        // Phase変更時に記事を再読み込み
        loadArticle();
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [events]);

  const loadArticle = async () => {
    try {
      setLoading(true);
      const data = await api.getArticle(articleId);
      setArticle(data);
      setError(null);
    } catch (err) {
      setError(err instanceof Error ? err.message : '記事の読み込みに失敗しました');
    } finally {
      setLoading(false);
    }
  };

  const handleExecutePhase = async (phase: string) => {
    setExecuting(phase);
    setError(null);
    try {
      let result: ArticleViewModel;
      switch (phase) {
        case 'plan':
          result = await api.executePlan(articleId);
          break;
        case 'do':
          result = await api.executeDo(articleId);
          break;
        case 'check':
          result = await api.executeCheck(articleId);
          break;
        case 'act':
          result = await api.executeAct(articleId);
          break;
        default:
          setExecuting(null);
          return;
      }
      setArticle(result);
      // 成功メッセージ
      const phaseNames: Record<string, string> = {
        plan: 'Plan Phase',
        do: 'Do Phase',
        check: 'Check Phase',
        act: 'Act Phase'
      };
      setSuccessMessage(`${phaseNames[phase] || phase}が正常に完了しました`);
      setTimeout(() => setSuccessMessage(null), 5000);
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'Phase実行に失敗しました';
      setError(errorMessage);
      console.error(`${phase} Phase実行エラー:`, err);
    } finally {
      setExecuting(null);
    }
  };

  const handlePublish = async () => {
    if (!confirm('Qiitaに投稿しますか？')) return;
    setError(null);
    setSuccessMessage(null);
    try {
      const result = await api.publishToQiita(articleId);
      await loadArticle();
      if (result.url) {
        setSuccessMessage(`投稿が完了しました: ${result.url}`);
      } else {
        setSuccessMessage('投稿が完了しました');
      }
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : '投稿に失敗しました';
      setError(errorMessage);
      console.error('投稿エラー:', err);
    }
  };

  const handleSaveContent = async (content: string) => {
    try {
      await api.updateArticle(articleId, { content });
      await loadArticle();
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : '保存に失敗しました';
      setError(errorMessage);
      console.error('保存エラー:', err);
    }
  };

  if (loading) {
    return (
      <div className="flex flex-col items-center justify-center min-h-[400px]">
        <LoadingSpinner />
        <p className="mt-4 text-gray-600 dark:text-gray-400">記事を読み込んでいます...</p>
      </div>
    );
  }

  if (error || !article) {
    return (
      <div className="p-4">
        <div className="p-4 bg-red-100 border border-red-400 text-red-700 rounded">
          <div className="flex justify-between items-center">
            <span>{error || '記事が見つかりません'}</span>
            <button
              onClick={() => router.push('/')}
              className="ml-4 px-3 py-1 bg-red-600 text-white rounded hover:bg-red-700 text-sm"
            >
              ダッシュボードに戻る
            </button>
          </div>
        </div>
      </div>
    );
  }

  const tabs = [
    {
      id: 'article',
      label: 'Article',
      content: (
        <ArticleViewer
          markdown={article.markdown || ''}
          editable={true}
          onSave={handleSaveContent}
        />
      ),
    },
    {
      id: 'process',
      label: 'AI Process',
      content: (
        <div>
          <h3 className="text-lg font-semibold mb-4">進捗タイムライン</h3>
          <div className="space-y-4">
            <div className="p-4 bg-gray-50 dark:bg-gray-800 rounded">
              <div className="flex justify-between items-center mb-2">
                <span className="font-medium">現在のPhase</span>
                <StatusBadge status={article.uiStatusText} phase={article.phase} />
              </div>
              <p className="text-sm text-gray-600 dark:text-gray-400">{article.nextActionHint}</p>
              {connected && (
                <div className="mt-2 text-xs text-green-600 dark:text-green-400">
                  ● リアルタイム更新中
                </div>
              )}
            </div>
            
            {/* ストリーミングイベント履歴 */}
            {events.length > 0 && (
              <div className="p-4 bg-gray-50 dark:bg-gray-800 rounded">
                <h4 className="text-sm font-medium mb-2">イベント履歴</h4>
                <div className="space-y-2 max-h-48 overflow-y-auto">
                  {events.slice(-10).reverse().map((event, idx) => (
                    <div key={idx} className="text-xs text-gray-600 dark:text-gray-400">
                      {event.type === 'phase_change' && (
                        <div>
                          Phase変更: {event.status_text} ({event.phase})
                        </div>
                      )}
                      {event.type === 'connected' && (
                        <div className="text-green-600 dark:text-green-400">
                          ストリーミング接続完了
                        </div>
                      )}
                      {event.type === 'error' && (
                        <div className="text-red-600 dark:text-red-400">
                          エラー: {event.message}
                        </div>
                      )}
                    </div>
                  ))}
                </div>
              </div>
            )}
            
            <div className="space-y-2">
              <button
                onClick={() => handleExecutePhase('plan')}
                disabled={executing !== null}
                className="w-full px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:opacity-50"
              >
                {executing === 'plan' ? '実行中...' : 'Plan Phase実行'}
              </button>
              <button
                onClick={() => handleExecutePhase('do')}
                disabled={executing !== null}
                className="w-full px-4 py-2 bg-yellow-600 text-white rounded-md hover:bg-yellow-700 disabled:opacity-50"
              >
                {executing === 'do' ? '実行中...' : 'Do Phase実行'}
              </button>
              <button
                onClick={() => handleExecutePhase('check')}
                disabled={executing !== null}
                className="w-full px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700 disabled:opacity-50"
              >
                {executing === 'check' ? '実行中...' : 'Check Phase実行'}
              </button>
              <button
                onClick={() => handleExecutePhase('act')}
                disabled={executing !== null}
                className="w-full px-4 py-2 bg-purple-600 text-white rounded-md hover:bg-purple-700 disabled:opacity-50"
              >
                {executing === 'act' ? '実行中...' : 'Act Phase実行'}
              </button>
            </div>
          </div>
        </div>
      ),
    },
    {
      id: 'feedback',
      label: 'Feedback',
      content: (
        <div className="space-y-6">
          <HumanFeedbackPanel articleId={articleId} onFeedbackAdded={loadArticle} />
          <div>
            <h3 className="text-lg font-semibold mb-4">フィードバック履歴</h3>
            <div className="space-y-4">
              {article.feedbackHistory && Array.isArray(article.feedbackHistory) && article.feedbackHistory.length > 0 ? (
                article.feedbackHistory.map((fb, idx) => (
                  <div key={idx} className="p-4 bg-gray-50 dark:bg-gray-800 rounded">
                    <p className="text-sm text-gray-600 dark:text-gray-400 mb-2">{fb.content || '（内容なし）'}</p>
                    <div className="text-xs text-gray-500">
                      {fb.created_at ? format(new Date(fb.created_at), 'yyyy年MM月dd日 HH:mm', { locale: ja }) : '-'}
                    </div>
                  </div>
                ))
              ) : (
                <p className="text-gray-600 dark:text-gray-400">フィードバックはありません</p>
              )}
            </div>
          </div>
        </div>
      ),
    },
    {
      id: 'analysis',
      label: 'Analysis',
      content: (
        <div className="space-y-6">
          <div>
            <div className="flex justify-between items-center mb-4">
              <h3 className="text-lg font-semibold">パフォーマンス分析</h3>
              <button
                onClick={() => handleExecutePhase('check')}
                disabled={executing !== null}
                className="px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700 disabled:opacity-50 text-sm"
              >
                {executing === 'check' ? '分析中...' : '分析を実行'}
              </button>
            </div>
            {article.analysisResults ? (
              <div className="p-4 bg-gray-50 dark:bg-gray-800 rounded">
                <h4 className="font-medium mb-2">分析結果</h4>
                <pre className="text-sm overflow-auto max-h-96">
                  {JSON.stringify(article.analysisResults, null, 2)}
                </pre>
              </div>
            ) : (
              <div className="p-4 bg-gray-50 dark:bg-gray-800 rounded">
                <p className="text-gray-600 dark:text-gray-400">分析結果はありません。分析を実行してください。</p>
              </div>
            )}
          </div>
          {article.kpiSummary && (
            <div className="p-4 bg-gray-50 dark:bg-gray-800 rounded">
              <h4 className="font-medium mb-2">KPIサマリ</h4>
              <div className="space-y-2">
                {typeof article.kpiSummary === 'object' ? (
                  <div className="space-y-2">
                    {Object.entries(article.kpiSummary).map(([key, value]) => (
                      <div key={key} className="flex justify-between items-center">
                        <span className="text-sm text-gray-600 dark:text-gray-400">{key}:</span>
                        <span className="text-sm font-medium">{String(value)}</span>
                      </div>
                    ))}
                  </div>
                ) : (
                  <pre className="text-sm overflow-auto max-h-96">
                    {JSON.stringify(article.kpiSummary, null, 2)}
                  </pre>
                )}
              </div>
            </div>
          )}
        </div>
      ),
    },
  ];

  return (
    <div>
      <div className="mb-8">
        <div className="flex justify-between items-start mb-4">
          <div>
            <h1 className="text-3xl font-bold text-gray-900 dark:text-white mb-2">{article.title}</h1>
            <StatusBadge status={article.uiStatusText} phase={article.phase} />
          </div>
          {article.qiitaUrl ? (
            <a
              href={article.qiitaUrl}
              target="_blank"
              rel="noopener noreferrer"
              className="px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700"
            >
              Qiitaで見る
            </a>
          ) : (
            <button
              onClick={handlePublish}
              className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
            >
              Qiitaに投稿
            </button>
          )}
        </div>
        <p className="text-gray-600 dark:text-gray-400 mb-2">{article.nextActionHint}</p>
        <div className="text-sm text-gray-500">
          {article.createdAt ? (
            <span>作成: {format(new Date(article.createdAt), 'yyyy年MM月dd日 HH:mm', { locale: ja })}</span>
          ) : null}
          {article.updatedAt ? (
            <span className="ml-4">
              更新: {format(new Date(article.updatedAt), 'yyyy年MM月dd日 HH:mm', { locale: ja })}
            </span>
          ) : null}
        </div>
      </div>

      {error && (
        <div className="mb-4 p-4 bg-red-100 border border-red-400 text-red-700 rounded">
          <div className="flex justify-between items-center">
            <span>{error}</span>
            <button
              onClick={() => setError(null)}
              className="text-red-700 hover:text-red-900"
              aria-label="エラーメッセージを閉じる"
            >
              ×
            </button>
          </div>
        </div>
      )}
      {successMessage && (
        <div className="mb-4 p-4 bg-green-100 border border-green-400 text-green-700 rounded">
          <div className="flex justify-between items-center">
            <span>{successMessage}</span>
            <button
              onClick={() => setSuccessMessage(null)}
              className="text-green-700 hover:text-green-900"
              aria-label="成功メッセージを閉じる"
            >
              ×
            </button>
          </div>
        </div>
      )}

      <div className="bg-white dark:bg-gray-800 rounded-lg shadow p-6">
        <TabContainer tabs={tabs} defaultTab="article" />
      </div>
    </div>
  );
}

