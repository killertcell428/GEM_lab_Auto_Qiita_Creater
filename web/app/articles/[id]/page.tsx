'use client';

import { useEffect, useState } from 'react';
import { useParams, useRouter } from 'next/navigation';
import { api, ArticleViewModel } from '@/lib/api';
import StatusBadge from '@/components/StatusBadge';
import ArticleViewer from '@/components/ArticleViewer';
import HumanFeedbackPanel from '@/components/HumanFeedbackPanel';
import LoadingSpinner from '@/components/LoadingSpinner';
import PerformanceChart from '@/components/PerformanceChart';
import WordCloud from '@/components/WordCloud';
import { useArticleStream } from '@/hooks/useArticleStream';
import { format } from 'date-fns';
import { ja } from 'date-fns/locale/ja';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

export default function ArticleDetailPage() {
  const params = useParams();
  const router = useRouter();
  const articleId = params.id as string;

  const [article, setArticle] = useState<ArticleViewModel | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [executing, setExecuting] = useState<string | null>(null);
  const [successMessage, setSuccessMessage] = useState<string | null>(null);
  const [keywords, setKeywords] = useState<Array<{ text: string; weight: number }>>([]);
  const [feedbackOpen, setFeedbackOpen] = useState<Record<string, boolean>>({
    plan: false,
    do: false,
    publish: false,
    check: false,
    act: false,
  });
  const scrollToSection = (id: string) => {
    const el = document.getElementById(id);
    if (el) el.scrollIntoView({ behavior: 'smooth', block: 'start' });
  };
  
  // ストリーミングフック
  const { events, connected } = useArticleStream(articleId);

  useEffect(() => {
    loadArticle();
    loadKeywords();
  }, [articleId]);
  
  const loadKeywords = async () => {
    try {
      const data = await api.getArticleKeywords(articleId);
      setKeywords(data.combined_keywords || []);
    } catch (err) {
      console.error('キーワード取得エラー:', err);
    }
  };
  
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

  // フェーズごとの簡易ステータスを算出
  const phaseOrder = ['plan', 'do', 'publish', 'check', 'act'];
  const phaseIndex = phaseOrder.indexOf(article?.phase || 'plan');
  const getPhaseStatus = (target: string) => {
    const idx = phaseOrder.indexOf(target);
    if (target === 'publish' && article?.pendingApproval) return '承認待ち';
    if (idx === -1) return '未実行';
    if (phaseIndex < idx) return '未実行';
    if (phaseIndex === idx) return '進行中';
    return '完了';
  };

  // Plan構成案をMarkdown風に整形
  const renderPlanMarkdown = (plan: any) => {
    if (!plan) return '';
    try {
      if (typeof plan === 'string') return plan;
      const title = plan.title || '構成案';
      const sections = Array.isArray(plan.sections) ? plan.sections : [];
      const lines = [`# ${title}`];
      sections.forEach((sec: any, i: number) => {
        const heading = sec.heading || `セクション${i + 1}`;
        lines.push(`## ${heading}`);
        if (sec.content_outline) lines.push(sec.content_outline);
        if (Array.isArray(sec.code_examples) && sec.code_examples.length > 0) {
          lines.push('### コード例');
          sec.code_examples.forEach((c: any) => lines.push(`- ${c}`));
        }
      });
      return lines.join('\n\n');
    } catch {
      return JSON.stringify(plan, null, 2);
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

  const handleApprove = async () => {
    if (!confirm('この記事を承認してQiitaに投稿しますか？')) return;
    setError(null);
    setSuccessMessage(null);
    try {
      const result = await api.approveArticle(articleId);
      await loadArticle();
      if (result.url) {
        setSuccessMessage(`承認が完了し、投稿しました: ${result.url}`);
      } else {
        setSuccessMessage('承認が完了しました');
      }
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : '承認に失敗しました';
      setError(errorMessage);
      console.error('承認エラー:', err);
    }
  };

  // 縦スクロールレイアウトのみ
  const renderFeedbackPanel = (phaseKey: string, defaultPhase: string) => (
    <div className="border border-dashed border-gray-200 dark:border-gray-700 rounded p-3 bg-gray-50 dark:bg-gray-900/50">
      <button
        type="button"
        onClick={() => setFeedbackOpen((prev) => ({ ...prev, [phaseKey]: !prev[phaseKey] }))}
        className="flex items-center justify-between w-full text-sm font-medium text-blue-700 dark:text-blue-300 hover:underline"
      >
        <span>＋ Human Feedback を追加</span>
        <span className="text-xs text-gray-500">{feedbackOpen[phaseKey] ? '閉じる' : '開く'}</span>
      </button>
      {feedbackOpen[phaseKey] && (
        <div className="mt-3">
          <HumanFeedbackPanel articleId={articleId} onFeedbackAdded={loadArticle} defaultPhase={defaultPhase} />
        </div>
      )}
    </div>
  );

  return (
    <div className="grid grid-cols-1 lg:grid-cols-[240px_1fr] gap-6">
      {/* サイドバー（固定） */}
      <aside className="hidden lg:block sticky top-4 self-start">
        <div className="bg-white dark:bg-gray-800 rounded-lg shadow p-4 space-y-3">
          <h3 className="text-sm font-semibold text-gray-700 dark:text-gray-200">フェーズナビ</h3>
          <div className="space-y-2">
            {[
              { id: 'overview', label: 'Overview', status: article.uiStatusText || '' },
              { id: 'plan', label: 'Plan', status: getPhaseStatus('plan') },
              { id: 'do', label: 'Do', status: getPhaseStatus('do') },
              { id: 'publish', label: 'Publish', status: article.pendingApproval ? '承認待ち' : getPhaseStatus('publish') },
              { id: 'check', label: 'Check', status: getPhaseStatus('check') },
              { id: 'act', label: 'Act', status: getPhaseStatus('act') },
              { id: 'feedback', label: 'Feedback', status: '履歴' },
            ].map((item) => (
              <button
                key={item.id}
                onClick={() => scrollToSection(item.id)}
                className="w-full flex items-center justify-between px-3 py-2 rounded hover:bg-gray-100 dark:hover:bg-gray-700 text-sm text-left"
              >
                <span>{item.label}</span>
                <span className="text-xs px-2 py-1 rounded bg-gray-100 dark:bg-gray-700 text-gray-700 dark:text-gray-200">
                  {item.status}
                </span>
              </button>
            ))}
          </div>
        </div>
      </aside>

      <div className="space-y-8">
        {/* ヘッダー */}
        <div className="mb-4">
          <div className="flex justify-between items-start mb-3">
            <div>
              <h1 className="text-3xl font-bold text-gray-900 dark:text-white mb-2">{article.title}</h1>
              <StatusBadge status={article.uiStatusText} phase={article.phase} />
            </div>
            {article.qiitaUrl ? (
              <a href={article.qiitaUrl} target="_blank" rel="noopener noreferrer" className="px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700">
                Qiitaで見る
              </a>
            ) : (
              <button onClick={handlePublish} className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700">
                Qiitaに投稿
              </button>
            )}
          </div>
          <p className="text-gray-600 dark:text-gray-400 mb-2">{article.nextActionHint}</p>
          <div className="text-sm text-gray-500">
            {article.createdAt && <span>作成: {format(new Date(article.createdAt), 'yyyy年MM月dd日 HH:mm', { locale: ja })}</span>}
            {article.updatedAt && <span className="ml-4">更新: {format(new Date(article.updatedAt), 'yyyy年MM月dd日 HH:mm', { locale: ja })}</span>}
          </div>
        </div>

        {error && (
          <div className="mb-4 p-4 bg-red-100 border border-red-400 text-red-700 rounded">
            <div className="flex justify-between items-center">
              <span>{error}</span>
              <button onClick={() => setError(null)} className="text-red-700 hover:text-red-900" aria-label="エラーメッセージを閉じる">
                ×
              </button>
            </div>
          </div>
        )}
        {successMessage && (
          <div className="mb-4 p-4 bg-green-100 border border-green-400 text-green-700 rounded">
            <div className="flex justify-between items-center">
              <span>{successMessage}</span>
              <button onClick={() => setSuccessMessage(null)} className="text-green-700 hover:text-green-900" aria-label="成功メッセージを閉じる">
                ×
              </button>
            </div>
          </div>
        )}

        {/* Overview */}
        <section className="bg-white dark:bg-gray-800 rounded-lg shadow p-6 space-y-4" id="overview">
          <h2 className="text-xl font-semibold">概要 / 次のアクション</h2>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <div className="p-4 bg-gray-50 dark:bg-gray-900 rounded space-y-2">
              <div className="flex items-center space-x-2">
                <span className="text-sm font-medium text-gray-600 dark:text-gray-400">ステータス</span>
                <StatusBadge status={article.uiStatusText} phase={article.phase} />
              </div>
              <div className="text-sm text-gray-700 dark:text-gray-200">次のアクション: {article.nextActionHint}</div>
              <div className="text-xs text-gray-500">作成: {article.createdAt ? format(new Date(article.createdAt), 'yyyy年MM月dd日 HH:mm', { locale: ja }) : '-'}</div>
              <div className="text-xs text-gray-500">更新: {article.updatedAt ? format(new Date(article.updatedAt), 'yyyy年MM月dd日 HH:mm', { locale: ja }) : '-'}</div>
              {article.qiitaUrl && (
                <a href={article.qiitaUrl} target="_blank" rel="noopener noreferrer" className="text-blue-600 dark:text-blue-400 hover:underline text-sm">
                  Qiitaで見る →
                </a>
              )}
            </div>
            <div className="p-4 bg-gray-50 dark:bg-gray-900 rounded space-y-2">
              <p className="text-sm font-semibold">人の操作ガイド</p>
              <ul className="text-xs text-gray-700 dark:text-gray-300 list-disc pl-4 space-y-1">
                <li>Plan/Doで記事を生成・編集</li>
                <li>Publishで承認して投稿（承認待ち時）</li>
                <li>投稿後7日でAnalyzeを自動実行、結果を確認</li>
                <li>フィードバックは各フェーズの入力欄で記載</li>
              </ul>
            </div>
          </div>
          {article.pendingApproval && article.approvalStatus === 'pending' && (
            <div className="p-4 bg-yellow-50 dark:bg-yellow-900/20 border border-yellow-200 dark:border-yellow-800 rounded">
              <div className="flex justify-between items-center mb-2">
                <span className="font-medium text-yellow-800 dark:text-yellow-200">承認待ち</span>
                <button onClick={handleApprove} className="px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700 text-sm">
                  承認して投稿
                </button>
              </div>
              {article.approvalDeadline && (
                <p className="text-sm text-yellow-700 dark:text-yellow-300">
                  承認期限: {new Date(article.approvalDeadline).toLocaleString('ja-JP')}
                </p>
              )}
            </div>
          )}
        </section>

        {/* Plan */}
        <section className="bg-white dark:bg-gray-800 rounded-lg shadow p-6 space-y-4" id="plan">
          <div className="flex justify-between items-center">
            <div className="space-y-1">
              <div className="flex items-center space-x-2">
                <h2 className="text-xl font-semibold">Plan（調査・構成）</h2>
                <span className="px-2 py-1 text-xs rounded bg-gray-100 dark:bg-gray-700 text-gray-700 dark:text-gray-200">
                  {getPhaseStatus('plan')}
                </span>
              </div>
              <p className="text-sm text-gray-600 dark:text-gray-400">AIが調査と構成案を作成します。必要ならフィードバックで要望を記載。</p>
            </div>
            <button onClick={() => handleExecutePhase('plan')} disabled={executing !== null} className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:opacity-50 text-sm">
              {executing === 'plan' ? '実行中...' : 'Planを実行'}
            </button>
          </div>
          <div className="p-4 bg-gray-50 dark:bg-gray-900 rounded space-y-3">
            <h3 className="text-sm font-medium">リサーチ結果</h3>
            {article.researchReport ? (
              <div className="prose dark:prose-invert max-w-none text-sm whitespace-pre-wrap">
                <ReactMarkdown remarkPlugins={[remarkGfm]}>{article.researchReport}</ReactMarkdown>
              </div>
            ) : (
              <p className="text-gray-600 dark:text-gray-400">まだリサーチ結果はありません。Planを実行してください。</p>
            )}
          </div>
          <div className="p-4 bg-gray-50 dark:bg-gray-900 rounded space-y-3">
            <h3 className="text-sm font-medium">構成案プレビュー</h3>
            {article.plan ? (
              <div className="prose dark:prose-invert max-w-none text-sm">
                <ReactMarkdown remarkPlugins={[remarkGfm]}>{renderPlanMarkdown(article.plan)}</ReactMarkdown>
              </div>
            ) : (
              <p className="text-gray-600 dark:text-gray-400">まだ構成案はありません。Planを実行してください。</p>
            )}
          </div>
          {renderFeedbackPanel('plan', 'plan')}
        </section>

        {/* Do */}
        <section className="bg-white dark:bg-gray-800 rounded-lg shadow p-6 space-y-4" id="do">
          <div className="flex justify-between items-center">
            <div className="space-y-1">
              <div className="flex items-center space-x-2">
                <h2 className="text-xl font-semibold">Do（執筆）</h2>
                <span className="px-2 py-1 text-xs rounded bg-gray-100 dark:bg-gray-700 text-gray-700 dark:text-gray-200">
                  {getPhaseStatus('do')}
                </span>
              </div>
              <p className="text-sm text-gray-600 dark:text-gray-400">Planの内容を基に本文を生成・編集します。</p>
            </div>
            <button onClick={() => handleExecutePhase('do')} disabled={executing !== null} className="px-4 py-2 bg-yellow-600 text-white rounded-md hover:bg-yellow-700 disabled:opacity-50 text-sm">
              {executing === 'do' ? '実行中...' : 'Doを実行'}
            </button>
          </div>
          <ArticleViewer markdown={article.markdown || ''} editable={true} onSave={handleSaveContent} />
          {renderFeedbackPanel('do', 'do')}
        </section>

        {/* Publish */}
        <section className="bg-white dark:bg-gray-800 rounded-lg shadow p-6 space-y-4" id="publish">
          <div className="flex justify-between items-center">
            <div className="space-y-1">
              <div className="flex items-center space-x-2">
                <h2 className="text-xl font-semibold">Publish（承認・投稿）</h2>
                <span className="px-2 py-1 text-xs rounded bg-gray-100 dark:bg-gray-700 text-gray-700 dark:text-gray-200">
                  {article.pendingApproval ? '承認待ち' : getPhaseStatus('publish')}
                </span>
              </div>
              <p className="text-sm text-gray-600 dark:text-gray-400">承認後にQiitaへ投稿します。</p>
            </div>
            <button onClick={handleApprove} disabled={executing !== null || !article.pendingApproval} className="px-4 py-2 bg-emerald-600 text-white rounded-md hover:bg-emerald-700 disabled:opacity-50 text-sm">
              承認して投稿
            </button>
          </div>
          {article.approvalDeadline && <p className="text-sm text-gray-500">承認期限: {new Date(article.approvalDeadline).toLocaleString('ja-JP')}</p>}
          {article.qiitaUrl ? (
            <div className="text-sm">
              <a href={article.qiitaUrl} target="_blank" rel="noopener noreferrer" className="text-blue-600 dark:text-blue-400 hover:underline">
                Qiitaで見る →
              </a>
            </div>
          ) : (
            <p className="text-gray-600 dark:text-gray-400 text-sm">承認待ちまたは未投稿です。</p>
          )}
          {renderFeedbackPanel('publish', 'publish')}
        </section>

        {/* Check */}
        <section className="bg-white dark:bg-gray-800 rounded-lg shadow p-6 space-y-4" id="check">
          <div className="flex justify-between items-center">
            <div className="space-y-1">
              <div className="flex items-center space-x-2">
                <h2 className="text-xl font-semibold">Check（分析）</h2>
                <span className="px-2 py-1 text-xs rounded bg-gray-100 dark:bg-gray-700 text-gray-700 dark:text-gray-200">
                  {getPhaseStatus('check')}
                </span>
              </div>
              <p className="text-sm text-gray-600 dark:text-gray-400">記事の分析を実行し、結果を確認します。</p>
            </div>
            <button onClick={() => handleExecutePhase('check')} disabled={executing !== null} className="px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700 disabled:opacity-50 text-sm">
              {executing === 'check' ? '実行中...' : 'Checkを実行'}
            </button>
          </div>
          <div className="p-4 bg-gray-50 dark:bg-gray-900 rounded">
            <PerformanceChart articleId={articleId} />
          </div>
          {article.analysisResults ? (
            <div className="p-4 bg-gray-50 dark:bg-gray-900 rounded">
              <h3 className="text-sm font-medium mb-2">分析結果</h3>
              <pre className="text-sm overflow-auto max-h-96">{JSON.stringify(article.analysisResults, null, 2)}</pre>
            </div>
          ) : (
            <p className="text-gray-600 dark:text-gray-400 text-sm">まだ分析結果はありません。投稿後、必要に応じて実行してください。</p>
          )}
          {article.kpiSummary && (
            <div className="p-4 bg-gray-50 dark:bg-gray-900 rounded">
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
                  <pre className="text-sm overflow-auto max-h-96">{JSON.stringify(article.kpiSummary, null, 2)}</pre>
                )}
              </div>
            </div>
          )}
          {keywords.length > 0 && (
            <div>
              <h3 className="text-lg font-semibold mb-4">キーワード分析</h3>
              <WordCloud words={keywords} />
            </div>
          )}
          {renderFeedbackPanel('check', 'check')}
        </section>

        {/* Act */}
        <section className="bg-white dark:bg-gray-800 rounded-lg shadow p-6 space-y-4" id="act">
          <div className="flex justify-between items-center">
            <div className="space-y-1">
              <div className="flex items-center space-x-2">
                <h2 className="text-xl font-semibold">Act（改善指針）</h2>
                <span className="px-2 py-1 text-xs rounded bg-gray-100 dark:bg-gray-700 text-gray-700 dark:text-gray-200">
                  {getPhaseStatus('act')}
                </span>
              </div>
              <p className="text-sm text-gray-600 dark:text-gray-400">分析結果を踏まえ、改善方針を決定します。</p>
            </div>
            <button onClick={() => handleExecutePhase('act')} disabled={executing !== null} className="px-4 py-2 bg-purple-600 text-white rounded-md hover:bg-purple-700 disabled:opacity-50 text-sm">
              {executing === 'act' ? '実行中...' : 'Actを実行'}
            </button>
          </div>
          {article.analysisResults ? (
            <div className="p-4 bg-gray-50 dark:bg-gray-900 rounded">
              <h3 className="text-sm font-medium mb-2">改善ヒント</h3>
              <p className="text-sm text-gray-700 dark:text-gray-200">分析結果をもとに、次のPlanに活かしてください。</p>
            </div>
          ) : (
            <p className="text-gray-600 dark:text-gray-400 text-sm">まだ分析が未実行のため、改善指針は表示できません。</p>
          )}
          {renderFeedbackPanel('act', 'act')}
        </section>

        {/* Feedbackまとめ */}
        <section className="bg-white dark:bg-gray-800 rounded-lg shadow p-6 space-y-4" id="feedback">
          <h2 className="text-xl font-semibold">Feedbackまとめ</h2>
          <div className="space-y-4">
            {article.feedbackHistory && Array.isArray(article.feedbackHistory) && article.feedbackHistory.length > 0 ? (
              article.feedbackHistory.map((fb, idx) => (
                <div key={idx} className="p-4 bg-gray-50 dark:bg-gray-900 rounded">
                  <p className="text-sm text-gray-600 dark:text-gray-400 mb-2">{fb.content || '（内容なし）'}</p>
                  <div className="text-xs text-gray-500">
                    {fb.created_at ? format(new Date(fb.created_at), 'yyyy年MM月dd日 HH:mm', { locale: ja }) : '-'}
                    {fb.phase && <span className="ml-2 px-2 py-1 bg-blue-100 text-blue-800 dark:bg-blue-900 dark:text-blue-200 rounded">{fb.phase}</span>}
                  </div>
                </div>
              ))
            ) : (
              <p className="text-gray-600 dark:text-gray-400">フィードバックはありません</p>
            )}
          </div>
        </section>
      </div>
    </div>
  );
}

